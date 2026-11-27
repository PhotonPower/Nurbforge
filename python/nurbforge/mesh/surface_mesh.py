"""Closed, consistently oriented surface meshes for boundary element solvers.

:meth:`GmshSession.surface_mesh` collects the triangles of all meshed surfaces into one
:class:`SurfaceMesh`, checks that the mesh is watertight (every edge shared by exactly two
triangles) and orients it: gmsh orients each surface by its parametrisation, so neighbouring
surfaces of a B-rep often disagree. The triangles are made consistent across surfaces and every
closed component is turned so that its normals point out of the solid – away from the material
for an outer boundary and into the void for a cavity. The orientation is applied to the gmsh
model as well (``gmsh.model.mesh.reverse``), so ``GmshSession.write`` exports the same mesh.
"""
from __future__ import annotations

from collections import deque
from dataclasses import dataclass, field

import numpy as np

# gmsh element types: 3-node and 6-node triangles
_TRIANGLE_TYPES = {2: 1, 9: 2}
# node permutation that reverses a triangle: (0, 2, 1) plus the matching edge midpoints
_REVERSE = {1: [0, 2, 1], 2: [0, 2, 1, 5, 4, 3]}


class MeshTopologyError(ValueError):
    """The surface mesh is not a closed, orientable 2-manifold."""


@dataclass
class SurfaceMesh:
    """A triangle surface mesh with 0-based connectivity.

    ``nodes`` (N, 3) coordinates; ``triangles`` (M, 3) for linear or (M, 6) for quadratic
    elements (gmsh node order: corners, then the midpoints of edges 0-1, 1-2, 2-0);
    ``surfaces`` (M,) the gmsh surface tag of every triangle; ``node_tags`` (N,) the gmsh node
    tags; ``groups`` triangle indices per face label (shells added with ``add_shell``).
    Meshes not made by gmsh (``shrink_wrap``) have the surface tag 0 and node tags 1..N."""

    nodes: np.ndarray
    triangles: np.ndarray
    surfaces: np.ndarray
    node_tags: np.ndarray
    groups: dict[str, np.ndarray] = field(default_factory=dict)

    @property
    def order(self) -> int:
        return 1 if self.triangles.shape[1] == 3 else 2

    def _corners(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        c = self.triangles[:, :3]
        return self.nodes[c[:, 0]], self.nodes[c[:, 1]], self.nodes[c[:, 2]]

    def normals(self) -> np.ndarray:
        """Unit normals of the (flat) corner triangles, (M, 3)."""
        a, b, c = self._corners()
        n = np.cross(b - a, c - a)
        return n / np.linalg.norm(n, axis=1, keepdims=True)

    def areas(self) -> np.ndarray:
        """Areas of the (flat) corner triangles, (M,)."""
        a, b, c = self._corners()
        return 0.5 * np.linalg.norm(np.cross(b - a, c - a), axis=1)

    def quality(self) -> np.ndarray:
        """Shape quality of the (flat) corner triangles, (M,): 4 sqrt(3) area / (sum of the
        squared side lengths), which is 1 for an equilateral triangle and 0 for a degenerate
        one."""
        a, b, c = self._corners()
        sides = ((b - a) ** 2).sum(axis=1) + ((c - b) ** 2).sum(axis=1) \
            + ((a - c) ** 2).sum(axis=1)
        return 4.0 * np.sqrt(3.0) * self.areas() / np.where(sides > 0.0, sides, 1.0)

    def centroids(self) -> np.ndarray:
        a, b, c = self._corners()
        return (a + b + c) / 3.0

    def volume(self) -> float:
        """Enclosed volume (divergence theorem over the corner triangles); positive for an
        outward-oriented closed mesh, cavities count negative."""
        a, b, c = self._corners()
        return float(np.einsum("ij,ij->i", a, np.cross(b, c)).sum() / 6.0)

    def boundary_edges(self) -> np.ndarray:
        """Corner-node edges used by only one triangle, (K, 2); empty for a watertight mesh."""
        edges, counts = _edge_counts(self.triangles)
        return edges[counts == 1]

    def is_watertight(self) -> bool:
        """Every edge is shared by exactly two triangles."""
        _, counts = _edge_counts(self.triangles)
        return bool(np.all(counts == 2))

    def is_consistently_oriented(self) -> bool:
        """Every directed edge occurs once, i.e. neighbours traverse a shared edge in opposite
        directions."""
        directed = _directed_edges(self.triangles)
        return len(np.unique(directed, axis=0)) == len(directed)

    def write(self, path) -> None:
        """Write the mesh to a file without gmsh: ``.msh`` (MSH 2.2 ASCII, with the surface
        tags as physical and elementary tags) or ``.stl`` (binary, corner triangles)."""
        path = str(path)
        ext = path.rsplit(".", 1)[-1].lower() if "." in path else ""
        if ext == "msh":
            etype = 2 if self.order == 1 else 9
            with open(path, "w", encoding="ascii") as f:
                f.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n")
                f.write(f"$Nodes\n{len(self.nodes)}\n")
                f.writelines(f"{tag} {x!r} {y!r} {z!r}\n" for tag, (x, y, z)
                             in zip(self.node_tags.tolist(), self.nodes.tolist()))
                f.write(f"$EndNodes\n$Elements\n{len(self.triangles)}\n")
                conn = self.node_tags[self.triangles].tolist()
                f.writelines(
                    f"{i} {etype} 2 {surface} {surface} {' '.join(map(str, nodes))}\n"
                    for i, (surface, nodes) in enumerate(zip(self.surfaces.tolist(), conn), 1))
                f.write("$EndElements\n")
        elif ext == "stl":
            a, b, c = self._corners()
            normal = np.cross(b - a, c - a)
            length = np.linalg.norm(normal, axis=1, keepdims=True)
            normal = normal / np.where(length > 0.0, length, 1.0)
            record = np.zeros(len(a), dtype=[("data", "<f4", (4, 3)), ("attribute", "<u2")])
            record["data"] = np.stack([normal, a, b, c], axis=1)
            with open(path, "wb") as f:
                f.write(b"Nurbforge surface mesh".ljust(80, b" "))
                f.write(np.uint32(len(a)).tobytes())
                f.write(record.tobytes())
        else:
            raise ValueError(f"SurfaceMesh.write: unsupported file type '{path}' "
                             "(use .msh or .stl)")


def _directed_edges(triangles: np.ndarray) -> np.ndarray:
    c = triangles[:, :3]
    return np.concatenate([c[:, [0, 1]], c[:, [1, 2]], c[:, [2, 0]]])


def _edge_counts(triangles: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    return np.unique(np.sort(_directed_edges(triangles), axis=1), axis=0, return_counts=True)


def _collect(gmsh) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Triangles of all meshed surfaces: (nodes, triangles, surface per triangle, node tags)."""
    tris, surf, order = [], [], None
    for _, tag in gmsh.model.getEntities(2):
        types, _, node_tags = gmsh.model.mesh.getElements(2, tag)
        for t, nt in zip(types, node_tags):
            if t not in _TRIANGLE_TYPES:
                raise MeshTopologyError(
                    f"surface {tag} has non-triangle elements (gmsh type {t}); "
                    "surface_mesh needs triangles (no recombination)")
            if order is not None and _TRIANGLE_TYPES[t] != order:
                raise MeshTopologyError("surface_mesh: mixed element orders")
            order = _TRIANGLE_TYPES[t]
            k = 3 * order
            block = np.asarray(nt, dtype=np.int64).reshape(-1, k)
            tris.append(block)
            surf.append(np.full(len(block), tag, dtype=np.int64))
    if not tris:
        raise MeshTopologyError("surface_mesh: no surface mesh – call mesh(dim=2) first")
    tri_tags = np.concatenate(tris)
    used = np.unique(tri_tags)
    all_tags, coords, _ = gmsh.model.mesh.getNodes()
    coords = np.asarray(coords).reshape(-1, 3)
    pos = {int(t): i for i, t in enumerate(all_tags)}
    nodes = coords[[pos[int(t)] for t in used]]
    triangles = np.searchsorted(used, tri_tags)
    return nodes, triangles, np.concatenate(surf), used


def _check_manifold(triangles: np.ndarray) -> None:
    _, counts = _edge_counts(triangles)
    open_edges, nonmanifold = int((counts == 1).sum()), int((counts > 2).sum())
    if open_edges or nonmanifold:
        parts = []
        if open_edges:
            parts.append(f"{open_edges} open edge(s) (used by one triangle)")
        if nonmanifold:
            parts.append(f"{nonmanifold} non-manifold edge(s) (used by more than two triangles)")
        raise MeshTopologyError(
            "surface mesh is not watertight: " + ", ".join(parts) + ". Mesh a closed shell "
            "added with add_shell, so that its faces share their edges.")


def _surface_flips(triangles: np.ndarray, surfaces: np.ndarray) -> tuple[dict[int, bool],
                                                                          list[list[int]]]:
    """Relative orientation of the surfaces: a flip per surface that makes all shared edges
    opposite, and the connected components (lists of surface tags)."""
    # each undirected edge is used by exactly two triangles (checked before); compare the
    # directions in which the two triangles traverse it
    c = triangles[:, :3]
    directed = np.concatenate([c[:, [0, 1]], c[:, [1, 2]], c[:, [2, 0]]])
    owner = np.tile(np.arange(len(c)), 3)
    key = np.sort(directed, axis=1)
    order = np.lexsort((key[:, 1], key[:, 0]))
    first, second = order[0::2], order[1::2]
    same_dir = directed[first, 0] == directed[second, 0]  # same direction: orientations differ
    sa, sb = surfaces[owner[first]], surfaces[owner[second]]

    # within a surface gmsh's triangles must already agree
    inner = sa == sb
    if np.any(same_dir & inner):
        raise MeshTopologyError("triangles of one gmsh surface are inconsistently oriented")
    # relation between surfaces: True = the two surfaces need opposite flips
    rel: dict[int, dict[int, bool]] = {}
    for a, b, s in zip(sa[~inner].tolist(), sb[~inner].tolist(), same_dir[~inner].tolist()):
        prev = rel.setdefault(a, {}).get(b)
        if prev is not None and prev != s:
            raise MeshTopologyError(
                f"surfaces {a} and {b} cannot be oriented consistently (non-orientable mesh)")
        rel[a][b] = s
        rel.setdefault(b, {})[a] = s
    flips: dict[int, bool] = {}
    components: list[list[int]] = []
    for start in dict.fromkeys(surfaces.tolist()):
        if start in flips:
            continue
        flips[start] = False
        comp, queue = [start], deque([start])
        while queue:
            s = queue.popleft()
            for t, opposite in rel.get(s, {}).items():
                want = flips[s] ^ opposite
                if t not in flips:
                    flips[t] = want
                    comp.append(t)
                    queue.append(t)
                elif flips[t] != want:
                    raise MeshTopologyError(
                        "surface mesh is not orientable (e.g. a Möbius-like face arrangement)")
        components.append(comp)
    return flips, components


def _signed_volume(nodes: np.ndarray, tri: np.ndarray) -> float:
    a, b, c = nodes[tri[:, 0]], nodes[tri[:, 1]], nodes[tri[:, 2]]
    return float(np.einsum("ij,ij->i", a, np.cross(b, c)).sum() / 6.0)


def _inside(point: np.ndarray, nodes: np.ndarray, tri: np.ndarray) -> bool:
    """Ray parity test: is ``point`` inside the closed triangle mesh ``tri``?"""
    d = np.array([0.5773, 0.6123, 0.5401])  # oblique direction, avoids axis-aligned edges
    a, b, c = nodes[tri[:, 0]], nodes[tri[:, 1]], nodes[tri[:, 2]]
    e1, e2 = b - a, c - a
    p = np.cross(d, e2)
    det = np.einsum("ij,ij->i", e1, p)
    ok = np.abs(det) > 1e-300
    inv = np.where(ok, 1.0 / np.where(ok, det, 1.0), 0.0)
    s = point - a
    u = np.einsum("ij,ij->i", s, p) * inv
    q = np.cross(s, e1)
    v = (q @ d) * inv
    t = np.einsum("ij,ij->i", e2, q) * inv
    hit = ok & (u >= 0) & (v >= 0) & (u + v <= 1) & (t > 0)
    return bool(hit.sum() % 2)


def build_surface_mesh(gmsh, labels: dict[str, list[int]] | None = None,
                       orient: bool = True) -> SurfaceMesh:
    """Collect, check and (with ``orient``) orient the surface mesh of the current gmsh model."""
    nodes, triangles, surfaces, node_tags = _collect(gmsh)
    _check_manifold(triangles)
    flips, components = _surface_flips(triangles, surfaces)
    if orient:
        rev = _REVERSE[triangles.shape[1] // 3]
        flip_tri = np.array([flips[s] for s in surfaces.tolist()])
        triangles[flip_tri] = triangles[flip_tri][:, rev]
        # outward: positive volume for a component nested in an even number of others
        comp_tris = [np.isin(surfaces, comp) for comp in components]
        to_reverse = set()
        for i, (comp, mask) in enumerate(zip(components, comp_tris)):
            corners = triangles[mask, :3]
            probe = nodes[corners[0]].mean(axis=0)
            depth = sum(_inside(probe, nodes, triangles[m, :3])
                        for j, m in enumerate(comp_tris) if j != i)
            vol = _signed_volume(nodes, corners)
            if (vol < 0) == (depth % 2 == 0):
                triangles[mask] = triangles[mask][:, rev]
                to_reverse.update(s for s in comp if not flips[s])
            else:
                to_reverse.update(s for s in comp if flips[s])
        if to_reverse:
            gmsh.model.mesh.reverse([(2, s) for s in sorted(to_reverse)])
    groups = {name: np.flatnonzero(np.isin(surfaces, tags))
              for name, tags in (labels or {}).items()}
    return SurfaceMesh(nodes, triangles, surfaces, node_tags, groups)
