"""Shrink wrap: one closed triangle mesh around arbitrary input (approximate, mesh level).

:func:`shrink_wrap` builds a watertight envelope around shells and loose surfaces, whether they
overlap, leave gaps or have naked edges. Unlike everything else in Nurbforge the result is
**not exact geometry**: it is a triangle mesh accurate to about the ``resolution``.

Method (numpy only): the input is sampled into points, the unsigned distance to them is
computed on a regular grid in a narrow band, the grid nodes that can be reached from outside
without coming closer to the input than the envelope distance are marked as outside, and the
surface between outside and the rest is extracted with marching tetrahedra, which always gives
a closed, consistently oriented 2-manifold.
"""
from __future__ import annotations

import itertools
import math
from collections.abc import Iterable

import numpy as np

from .._core import NurbsSurface, Shell, TrimmedSurface
from .surface_mesh import SurfaceMesh

# Samples lie at most _SPACING * resolution apart (along the grid lines of a face).
_SPACING = 0.25
# Smallest envelope distance that separates inside from outside on the grid, in units of the
# resolution: a grid edge crossing the input has an end within 0.5 of the crossing, which lies
# within _SPACING of a sample.
_WALL = 0.75
# Triangles of the envelope are sampled with this many intervals per side.
_SUBDIVISION = 3
# A part of the result must reach this deep below its surface (in units of the resolution).
_CORE = 0.1
# Step of the smoothing, as a fraction of the way to the centre of the neighbours.
_RELAXATION = 0.6


def _spacing(points: np.ndarray, axis: int) -> float:
    if points.shape[axis] < 2:
        return 0.0
    return float(np.linalg.norm(np.diff(points, axis=axis), axis=-1).max())


def _sample_surface(surface: NurbsSurface, step: float) -> tuple[np.ndarray, int, int]:
    """Grid of points on the surface with neighbours at most ``step`` apart."""
    nu = nv = 9
    for _ in range(8):
        pts = np.asarray(surface.sample(nu, nv))
        du, dv = _spacing(pts, 0), _spacing(pts, 1)
        if du <= step and dv <= step:
            break
        # the chords of a coarse grid underestimate the spacing of a fine one: add a margin
        if du > step:
            nu = math.ceil((nu - 1) * du / step * 1.2) + 1
        if dv > step:
            nv = math.ceil((nv - 1) * dv / step * 1.2) + 1
        if nu * nv > 50_000_000:
            raise ValueError("shrink_wrap: resolution too fine for the size of the input")
    return pts, nu, nv


def _sample_curve(curve, step: float) -> np.ndarray:
    n = 9
    for _ in range(8):
        pts = np.asarray(curve.sample(n))
        d = _spacing(pts, 0)
        if d <= step:
            break
        n = math.ceil((n - 1) * d / step * 1.2) + 1
    return pts


def _sample_face(face: TrimmedSurface, step: float) -> list[np.ndarray]:
    surface = face.surface
    pts, nu, nv = _sample_surface(surface, step)
    if face.is_untrimmed:
        return [pts.reshape(-1, 3)]
    (u0, u1), (v0, v1) = surface.domain_u, surface.domain_v
    us, vs = np.linspace(u0, u1, nu), np.linspace(v0, v1, nv)
    contains = face.contains
    keep = np.fromiter((contains(u, v) for u in us.tolist() for v in vs.tolist()),
                       dtype=bool, count=nu * nv)
    out = [pts.reshape(-1, 3)[keep]]
    # the grid stops up to one step short of the trim curves
    for loop in face.boundary_curves():
        out.extend(_sample_curve(c, step) for c in loop)
    return out


def _sample(objects, step: float) -> np.ndarray:
    if isinstance(objects, (Shell, TrimmedSurface, NurbsSurface)):
        objects = [objects]
    elif not isinstance(objects, Iterable):
        raise TypeError("shrink_wrap: expected shells or surfaces")
    clouds: list[np.ndarray] = []
    for obj in objects:
        if isinstance(obj, Shell):
            for face in obj.faces:
                clouds.extend(_sample_face(face, step))
        elif isinstance(obj, TrimmedSurface):
            clouds.extend(_sample_face(obj, step))
        elif isinstance(obj, NurbsSurface):
            clouds.append(_sample_surface(obj, step)[0].reshape(-1, 3))
        else:
            raise TypeError(f"shrink_wrap: cannot wrap a {type(obj).__name__}; "
                            "expected Shell, TrimmedSurface or NurbsSurface")
    clouds = [c for c in clouds if len(c)]
    if not clouds:
        raise ValueError("shrink_wrap: nothing to wrap")
    return np.concatenate(clouds)


class _Grid:
    """Regular grid with one layer of padding; nodes are addressed by their flat index."""

    def __init__(self, lo: np.ndarray, hi: np.ndarray, h: float, margin: float,
                 max_nodes: int) -> None:
        self.h = h
        # +1 node of padding on every side, never touched by the algorithms below
        self.origin = lo - margin - h
        counts = np.ceil((hi - lo + 2.0 * margin) / h).astype(np.int64) + 3
        if int(np.prod(counts)) > max_nodes:
            raise ValueError(
                f"shrink_wrap: the grid would have {int(np.prod(counts)):,} nodes "
                f"({counts[0]} x {counts[1]} x {counts[2]}), more than max_nodes = "
                f"{max_nodes:,}; use a coarser resolution or raise max_nodes")
        self.shape = (int(counts[0]), int(counts[1]), int(counts[2]))
        self.size = self.shape[0] * self.shape[1] * self.shape[2]
        self.strides = np.array([self.shape[1] * self.shape[2], self.shape[2], 1],
                                dtype=np.int64)
        pad = np.ones(self.shape, dtype=bool)
        pad[1:-1, 1:-1, 1:-1] = False
        self.padding = pad.ravel()

    def position(self, index: np.ndarray) -> np.ndarray:
        ijk = np.stack(np.unravel_index(index, self.shape), axis=-1)
        return self.origin + self.h * ijk


def _distance_field(grid: _Grid, points: np.ndarray, band: float
                    ) -> tuple[np.ndarray, np.ndarray]:
    """Squared distance to the nearest sample and its index for all nodes within ``band`` of
    the samples (inf / -1 elsewhere), by propagating nearest samples from node to node."""
    d2 = np.full(grid.size, np.inf)
    near = np.full(grid.size, -1, dtype=np.int64)
    # seeds: the eight nodes around every sample
    cell = np.floor((points - grid.origin) / grid.h).astype(np.int64)
    rest = points - (grid.origin + grid.h * cell)
    base = cell @ grid.strides
    for corner in itertools.product((0, 1), repeat=3):
        node = base + int(np.dot(corner, grid.strides))
        dist = ((rest - grid.h * np.array(corner)) ** 2).sum(axis=-1)
        np.minimum.at(d2, node, dist)
        hit = dist == d2[node]
        near[node[hit]] = np.flatnonzero(hit)
    changed = np.flatnonzero(near >= 0)
    # a later corner may have lowered the distance of a node below that of its sample
    d2[changed] = ((grid.position(changed) - points[near[changed]]) ** 2).sum(axis=-1)

    shifts = np.array([s for s in itertools.product((-1, 0, 1), repeat=3) if any(s)],
                      dtype=np.int64) @ grid.strides
    limit = band * band
    blocked = grid.padding
    while len(changed):
        improved = []
        for shift in shifts.tolist():
            target = changed + shift
            cand = near[changed]
            ok = ~blocked[target]
            target, cand = target[ok], cand[ok]
            dist = ((grid.position(target) - points[cand]) ** 2).sum(axis=-1)
            ok = (dist < d2[target]) & (dist <= limit)
            target = target[ok]
            d2[target] = dist[ok]
            near[target] = cand[ok]
            improved.append(target)
        changed = np.unique(np.concatenate(improved))
    return d2, near


def _flood(grid: _Grid, free: np.ndarray, start: np.ndarray | None = None) -> np.ndarray:
    """Free nodes connected to the ``start`` nodes (default: the border of the grid) through
    free nodes (6-neighbourhood)."""
    free = free & ~grid.padding
    if start is None:
        start = np.array([int(grid.strides.sum())], dtype=np.int64)  # node (1, 1, 1)
    frontier = start[free[start]]
    reached = np.zeros(grid.size, dtype=bool)
    free[frontier] = False
    reached[frontier] = True
    shifts = np.concatenate([grid.strides, -grid.strides])
    while len(frontier):
        nb = (frontier[:, None] + shifts[None, :]).ravel()
        nb = np.unique(nb[free[nb]])
        free[nb] = False
        reached[nb] = True
        frontier = nb
    return reached


# Six tetrahedra per cube around the diagonal (0,0,0)-(1,1,1): the same split in every cube,
# so the faces of neighbouring cubes match.
_TETS = np.array([[(0, 0, 0), tuple(np.eye(3, dtype=int)[a]),
                   tuple(np.eye(3, dtype=int)[a] + np.eye(3, dtype=int)[b]), (1, 1, 1)]
                  for a, b, _ in itertools.permutations(range(3))], dtype=np.int64)
_OTHERS = np.array([[1, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]])


def _marching_tetrahedra(grid: _Grid, outside: np.ndarray
                         ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Surface between the outside nodes and the others.

    Returns the grid edges that carry a vertex as pairs (outside node, other node), the
    triangles as indices into them and, per triangle, a direction pointing outside."""
    inside = ~outside.reshape(grid.shape)
    count = np.zeros(tuple(n - 1 for n in grid.shape), dtype=np.uint8)
    for di, dj, dk in itertools.product((0, 1), repeat=3):
        count += inside[di:grid.shape[0] - 1 + di, dj:grid.shape[1] - 1 + dj,
                        dk:grid.shape[2] - 1 + dk]
    active = np.argwhere((count > 0) & (count < 8)) @ grid.strides
    inside = inside.ravel()

    pairs: list[np.ndarray] = []  # (T, 3, 2) grid nodes of the edges carrying the vertices
    towards: list[np.ndarray] = []  # (T, 3) a direction from inside to outside
    for tet in _TETS:
        nodes = active[:, None] + (tet @ grid.strides)[None, :]
        ins = inside[nodes]
        k = ins.sum(axis=1)
        for lone_inside in (True, False):
            rows = np.flatnonzero(k == (1 if lone_inside else 3))
            if not len(rows):
                continue
            n, i = nodes[rows], ins[rows]
            lone = np.argmax(i == lone_inside, axis=1)
            a = np.take_along_axis(n, lone[:, None], axis=1)
            b = np.take_along_axis(n, _OTHERS[lone], axis=1)
            pairs.append(np.stack([np.broadcast_to(a, b.shape), b], axis=-1))
            g = grid.position(b).mean(axis=1) - grid.position(a[:, 0])
            towards.append(g if lone_inside else -g)
        rows = np.flatnonzero(k == 2)
        if len(rows):
            n, i = nodes[rows], ins[rows]
            order = np.argsort(~i, axis=1, kind="stable")  # the two inside corners first
            n = np.take_along_axis(n, order, axis=1)
            a, b, c, d = n[:, 0], n[:, 1], n[:, 2], n[:, 3]
            ac, ad = np.stack([a, c], -1), np.stack([a, d], -1)
            bd, bc = np.stack([b, d], -1), np.stack([b, c], -1)
            g = (grid.position(c) + grid.position(d) - grid.position(a) - grid.position(b))
            for tri in ((ac, ad, bd), (ac, bd, bc)):
                pairs.append(np.stack(tri, axis=1))
                towards.append(g)
    if not pairs:
        raise ValueError("shrink_wrap: the envelope is empty")
    pair = np.concatenate(pairs)
    g = np.concatenate(towards)

    pair = np.sort(pair, axis=-1)
    key = pair[..., 0] * grid.size + pair[..., 1]
    _, first, tri = np.unique(key.ravel(), return_index=True, return_inverse=True)
    ends = pair.reshape(-1, 2)[first]
    swap = ~outside[ends[:, 0]]
    ends[swap] = ends[swap][:, ::-1]
    return ends, tri.reshape(-1, 3), g


def _place(grid: _Grid, ends: np.ndarray, t: np.ndarray, tri: np.ndarray, towards: np.ndarray
           ) -> tuple[np.ndarray, np.ndarray]:
    """Vertices at the fraction ``t`` of the way from the outside node of each edge to the
    other one, and the triangles turned so that their normals point outside."""
    t = np.clip(t, 1e-3, 1.0 - 1e-3)[:, None]
    pa, pb = grid.position(ends[:, 0]), grid.position(ends[:, 1])
    vertices = pa + t * (pb - pa)
    p = vertices[tri]
    normal = np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0])
    flip = np.einsum("ij,ij->i", normal, towards) < 0
    tri = tri.copy()
    tri[flip] = tri[flip][:, ::-1]
    return vertices, tri


def _relax(mesh: SurfaceMesh, min_quality: float, max_iterations: int = 50) -> None:
    """Move the nodes within the surface towards the centres of their neighbours until the
    worst triangle reaches ``min_quality`` (tangential Laplacian smoothing, in place). The
    connectivity is kept, so the mesh stays closed; the shape changes to second order only,
    since the nodes move in their tangent planes."""
    tri = mesh.triangles
    count = len(mesh.nodes)
    # closed mesh: every edge occurs once per direction, so these are all neighbours
    edge = np.concatenate([tri[:, [0, 1]], tri[:, [1, 2]], tri[:, [2, 0]]])
    degree = np.bincount(edge[:, 0], minlength=count).astype(float)[:, None]
    for iteration in range(max_iterations):
        if iteration % 5 == 0 and mesh.quality().min() >= min_quality:
            break
        v = mesh.nodes
        a, b, c = v[tri[:, 0]], v[tri[:, 1]], v[tri[:, 2]]
        face = np.cross(b - a, c - a)  # area-weighted face normals
        normal = np.stack([np.bincount(tri.ravel(), weights=np.repeat(face[:, k], 3),
                                       minlength=count) for k in range(3)], axis=1)
        normal /= np.maximum(np.linalg.norm(normal, axis=1, keepdims=True), 1e-300)
        centre = np.stack([np.bincount(edge[:, 0], weights=v[edge[:, 1], k], minlength=count)
                           for k in range(3)], axis=1) / degree
        step = centre - v
        step -= normal * np.einsum("ij,ij->i", step, normal)[:, None]
        mesh.nodes = v + _RELAXATION * step


def shrink_wrap(objects, resolution: float, offset: float = 0.0, fill_holes: bool = True,
                max_nodes: int = 20_000_000, min_quality: float = 0.3) -> SurfaceMesh:
    """One closed triangle mesh around shells and surfaces (approximate).

    ``objects``: a ``Shell``, ``TrimmedSurface`` or ``NurbsSurface``, or a list of them. They
    may overlap, leave gaps or be open. ``resolution`` is the cell size of the grid the envelope
    is built on, which is about the edge length of the mesh and the accuracy of the result:
    openings and gaps narrower than about 1.5 x ``resolution`` (plus twice the ``offset``) are
    closed, smaller details are smoothed away. ``offset`` is the distance of the envelope from
    the input.

    The envelope is first built at a distance of at least 0.75 x ``resolution``, which is what
    makes it watertight on the grid; for a smaller ``offset`` the surface that lies the
    difference further in is taken (a morphological closing). Parts thinner than about 1.5 x
    ``resolution`` do not survive this: an open surface needs an ``offset`` of at least
    0.75 x ``resolution``, a ``ValueError`` is raised if nothing is left.

    With ``fill_holes`` (default) everything that cannot be reached from outside is filled,
    and the mesh is the outer envelope only; without it, enclosed voids larger than the
    envelope distance are kept as inner components (normals pointing into the void). Note that
    the inside of a closed shell is such a void: the input is a set of surfaces.

    The triangles cut out of the grid can be arbitrarily thin. If the worst of them has a
    :meth:`SurfaceMesh.quality` below ``min_quality``, the nodes are moved within the surface
    until it is reached (at most 50 smoothing steps; typically 5 give a worst quality above
    0.3 and a mean of 0.85). ``min_quality=0`` returns the mesh as cut out of the grid. For
    triangles of a different size see :func:`nurbforge.mesh.remesh`.

    Returns a :class:`SurfaceMesh` with linear triangles, watertight and with normals pointing
    out of the wrapped volume. Parts further apart than the gap size give separate closed
    components. Raises ``ValueError`` if the grid would exceed ``max_nodes`` nodes.
    """
    if not resolution > 0.0:
        raise ValueError("shrink_wrap: resolution must be positive")
    if offset < 0.0:
        raise ValueError("shrink_wrap: offset must not be negative")
    h = float(resolution)
    iso = max(float(offset), _WALL * h)
    band = iso + 2.0 * h
    coarse = _sample(objects, math.inf)  # fails early if the grid would be too large
    _Grid(coarse.min(axis=0), coarse.max(axis=0), h, band, max_nodes)
    points = _sample(objects, _SPACING * h)
    grid = _Grid(points.min(axis=0), points.max(axis=0), h, band, max_nodes)
    d2, near = _distance_field(grid, points, band)

    free = d2 > iso * iso
    outside = _flood(grid, free) if fill_holes else (free & ~grid.padding)
    outside |= grid.padding
    ends, triangles, towards = _marching_tetrahedra(grid, outside)
    # Where the distance falls to iso on the way from the outside node to the other one. The
    # distance is unsigned: if the two nodes see the input in opposite directions it lies
    # between them, and the distance of the second node counts negative.
    out, other = ends[:, 0], ends[:, 1]
    d_out = np.minimum(np.sqrt(d2[out]), band)
    d_other = np.minimum(np.sqrt(d2[other]), band)
    to_out = grid.position(out) - points[near[out]]
    to_other = grid.position(other) - points[near[other]]
    across = (np.einsum("ij,ij->i", to_out, to_other) < 0.0) & (near[out] >= 0) \
        & (near[other] >= 0)
    d_other = np.where(across, -d_other, d_other)
    span = d_out - d_other
    t = np.where(span > 0.0, (d_out - iso) / np.where(span > 0.0, span, 1.0), 0.5)
    vertices, triangles = _place(grid, ends, t, triangles, towards)

    if iso > offset:
        # Closing: take the surface at the distance iso - offset inside the envelope. Away from
        # gaps and concave corners that is the surface at the distance offset from the input.
        back = iso - offset
        corner = vertices[triangles]
        weights = np.array([(i, j, _SUBDIVISION - i - j) for i in range(_SUBDIVISION + 1)
                            for j in range(_SUBDIVISION + 1 - i)], dtype=float) / _SUBDIVISION
        skin = np.einsum("kc,tcx->tkx", weights, corner).reshape(-1, 3)
        band = back + 2.0 * h
        d2, _ = _distance_field(grid, skin, band)
        dist = np.minimum(np.sqrt(d2), band)
        value = np.where(outside, dist + back, back - dist)  # signed distance, positive outside
        # Keep what is connected to a node clearly below the surface: a part thinner than the
        # envelope leaves only noise around zero, which would give scattered fragments.
        solid = _flood(grid, value <= 0.0, np.flatnonzero(value < -_CORE * h))
        outside = ~solid
        value = np.where(solid, value, np.maximum(value, 1e-3 * h))
        if outside.all():
            raise ValueError(
                "shrink_wrap: nothing is left of the envelope; the input is thinner than the "
                "resolution (open surfaces need an offset of at least 0.75 x resolution)")
        ends, triangles, towards = _marching_tetrahedra(grid, outside)
        v_out, v_other = value[ends[:, 0]], value[ends[:, 1]]
        vertices, triangles = _place(grid, ends, v_out / (v_out - v_other), triangles, towards)

    hull = SurfaceMesh(vertices, triangles, np.zeros(len(triangles), dtype=np.int64),
                       np.arange(1, len(vertices) + 1, dtype=np.int64))
    if min_quality > 0.0:
        _relax(hull, min_quality)
    return hull
