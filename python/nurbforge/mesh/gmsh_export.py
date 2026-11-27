"""Transfer Nurbforge NURBS geometry into gmsh's OpenCASCADE kernel and mesh it.

The NURBS data is handed over *exactly* (degree, knots, multiplicities, weights),
so gmsh meshes the true geometry rather than a tessellation.

Example
-------
>>> import nurbforge as nf
>>> from nurbforge.mesh import GmshSession
>>> sphere = nf.revolve(nf.make_arc((0,0,0),(0,0,-1),(1,0,0),1.0,0.0,3.14159265), (0,0,0),(0,0,1),6.2831853)
>>> with GmshSession() as g:                         # doctest: +SKIP
...     tag = g.add_surface(sphere)
...     g.mesh(dim=2, size=0.1)
...     g.write("sphere.msh")
"""
from __future__ import annotations

import math
import os
import tempfile
from collections.abc import Iterable

import numpy as np

from .._core import NurbsCurve, NurbsSurface, Shell, TrimmedSurface, unique_knots, write_step
from .surface_mesh import SurfaceMesh, build_surface_mesh


def _gmsh():
    try:
        import gmsh  # type: ignore
    except ImportError as exc:  # pragma: no cover
        raise ImportError(
            "gmsh is not installed. Install with `pip install nurbforge[mesh]` "
            "or `pip install gmsh`."
        ) from exc
    return gmsh


def _remove_control_points(occ, tags: list[int]) -> None:
    """Control points are only needed to build the spline; left in the model they
    would become free vertices and end up as stray mesh nodes."""
    occ.remove([(0, t) for t in tags])


def add_curve(curve: NurbsCurve) -> int:
    """Add a NURBS curve to the current gmsh OCC model. Returns the curve tag."""
    gmsh = _gmsh()
    occ = gmsh.model.occ
    tags = [occ.addPoint(float(p[0]), float(p[1]), float(p[2])) for p in curve.control_points]
    knots, mults = unique_knots(curve.knots)
    tag = occ.addBSpline(
        tags, degree=curve.degree, weights=list(curve.weights), knots=knots, multiplicities=mults
    )
    _remove_control_points(occ, tags)
    return tag


def add_surface(surface: NurbsSurface, wire_tags: list[int] | None = None) -> int:
    """Add a NURBS surface to the current gmsh OCC model. Returns the surface tag.

    gmsh expects control points with the u index running fastest, while Nurbforge
    stores them row-major with v fastest, hence the reordering below. ``wire_tags``
    (outer contour first, then holes) trims the surface with wires defined in its
    parameter space.
    """
    gmsh = _gmsh()
    occ = gmsh.model.occ
    nu, nv = surface.num_u, surface.num_v
    pts = surface.control_points
    w = surface.weights
    tags, weights = [], []
    for j in range(nv):
        for i in range(nu):
            p = pts[i * nv + j]
            tags.append(occ.addPoint(float(p[0]), float(p[1]), float(p[2])))
            weights.append(w[i * nv + j])
    ku, mu = unique_knots(surface.knots_u)
    kv, mv = unique_knots(surface.knots_v)
    tag = occ.addBSplineSurface(
        tags, nu, degreeU=surface.degree_u, degreeV=surface.degree_v, weights=weights,
        knotsU=ku, knotsV=kv, multiplicitiesU=mu, multiplicitiesV=mv,
        wireTags=wire_tags or [], wire3D=False,
    )
    _remove_control_points(occ, tags)
    return tag


def add_trimmed_surface(trimmed: TrimmedSurface) -> int:
    """Add a trimmed NURBS surface: the trim loops are passed exactly as B-spline wires in
    the parameter space of the surface (outer contour first, then holes). Returns the tag."""
    if trimmed.is_untrimmed:
        # no wires: OCC handles the full (possibly periodic) B-spline patch itself, whereas a
        # parameter-space rectangle would contain a closed surface's seam twice
        return add_surface(trimmed.surface)
    gmsh = _gmsh()
    occ = gmsh.model.occ
    curve_tags: list[int] = []
    wires: list[int] = []
    for loop in [trimmed.outer, *trimmed.holes]:
        tags = [add_curve(c) for c in loop.curves]
        curve_tags += tags
        wires.append(occ.addWire(tags))
    tag = add_surface(trimmed.surface, wires)
    # the parameter-space curves are only construction helpers
    occ.remove([(1, t) for t in curve_tags], recursive=True)
    return tag


class GmshSession:
    """Context manager around ``gmsh.initialize()`` / ``gmsh.finalize()``.

    Shells added with :meth:`add_shell` keep their face labels, which :meth:`face_tags` maps to
    gmsh surfaces and :meth:`refine` uses for local mesh sizes."""

    def __init__(self, model_name: str = "nurbforge", verbose: bool = False):
        self.model_name = model_name
        self.verbose = verbose
        self.gmsh = _gmsh()
        # (label, interior point) of every face added with add_shell; matched to the gmsh surfaces
        # only when needed, because removeAllDuplicates renumbers surfaces
        self._faces: list[tuple[str, tuple[float, float, float]]] = []
        self._refinements: list[tuple[list[str | int], float, float]] = []
        self._fields: list[int] = []

    def __enter__(self) -> GmshSession:  # noqa: PYI034
        g = self.gmsh
        if not g.isInitialized():
            g.initialize()
        g.option.setNumber("General.Terminal", 1 if self.verbose else 0)
        g.model.add(self.model_name)
        return self

    def __exit__(self, *exc) -> None:
        self.gmsh.finalize()

    def add_curve(self, curve: NurbsCurve) -> int:
        return add_curve(curve)

    def add_surface(self, surface: NurbsSurface | TrimmedSurface) -> int:
        if isinstance(surface, TrimmedSurface):
            return add_trimmed_surface(surface)
        return add_surface(surface)

    def add_shell(self, shell: Shell, solid: bool = True, physical_groups: bool = False,
                  volume_name: str = "volume") -> list[int]:
        """As the module function :func:`add_shell`; in addition the session remembers the face
        labels for :meth:`face_tags` and :meth:`refine`."""
        tags = add_shell(shell, solid=solid, physical_groups=physical_groups,
                         volume_name=volume_name)
        self._faces += [(f.label or "unlabeled", _interior_point(f)) for f in shell.faces]
        return tags

    def face_tags(self) -> dict[str, list[int]]:
        """The current gmsh surface tags of every face label of the shells added with
        :meth:`add_shell` (unlabelled faces: ``"unlabeled"``)."""
        g = self.gmsh
        g.model.occ.synchronize()
        tags = [t for _, t in g.model.getEntities(2)]
        out: dict[str, set[int]] = {}
        for label, p in self._faces:
            out.setdefault(label, set()).add(_closest_surface(g, tags, p))
        return {k: sorted(v) for k, v in out.items()}

    def refine(self, where: str | int | Iterable[str | int], size: float,
               distance: float = 0.0) -> None:
        """Local mesh size ``size`` on the given faces: face labels of shells added with
        :meth:`add_shell` and/or gmsh surface tags. With ``distance > 0`` the size grows linearly
        from ``size`` on the faces to the global size (``mesh(size=...)``) at that distance, also
        into the volume; with ``distance = 0`` only the faces themselves are refined. Applied by
        :meth:`mesh`; several refinements combine by taking the smallest size."""
        if not size > 0:
            raise ValueError("refine: size must be positive")
        if not distance >= 0:
            raise ValueError("refine: distance must be non-negative")
        items = [where] if isinstance(where, (str, int)) else list(where)
        if not items:
            raise ValueError("refine: no faces given")
        self._refinements.append((items, float(size), float(distance)))

    def mesh(self, dim: int = 2, size: float | None = None, order: int = 1,
             heal: bool = True, curvature: int | None = None) -> None:
        """Synchronise the OCC model and generate a mesh. ``size`` is the global maximum element
        size, ``curvature`` the number of elements per 2*pi of curvature (smaller elements on
        strongly curved faces), local sizes come from :meth:`refine`."""
        g = self.gmsh
        if heal:
            # merges duplicate control-point vertices and fixes seams/poles
            g.model.occ.removeAllDuplicates()
        g.model.occ.synchronize()
        if size is not None:
            g.option.setNumber("Mesh.MeshSizeMax", size)
        if curvature is not None:
            if curvature < 1:
                raise ValueError("mesh: curvature must be at least 1 element per 2*pi")
            g.option.setNumber("Mesh.MeshSizeFromCurvature", int(curvature))
        self._apply_refinements(size)
        g.model.mesh.generate(dim)
        if order > 1:
            g.model.mesh.setOrder(order)

    def _apply_refinements(self, size: float | None) -> None:
        """(Re)build the background size field from the refinements: per refinement a Constant
        field on the faces (exact size there) and, with a distance, a Threshold field on the
        distance to them; the background field is their minimum."""
        field = self.gmsh.model.mesh.field
        for f in self._fields:
            field.remove(f)
        self._fields = []
        if not self._refinements:
            return
        labels = self.face_tags() if any(isinstance(w, str) for r in self._refinements
                                         for w in r[0]) else {}
        if size is None:  # growth target without a global size: the model's diagonal
            x0, y0, z0, x1, y1, z1 = self.gmsh.model.getBoundingBox(-1, -1)
            size = math.dist((x0, y0, z0), (x1, y1, z1))
        mins = []
        for items, local, distance in self._refinements:
            surfaces: set[int] = set()
            for w in items:
                if isinstance(w, str):
                    if w not in labels:
                        raise KeyError(f"refine: no face labelled {w!r} "
                                       f"(labels: {', '.join(sorted(labels)) or 'none'})")
                    surfaces.update(labels[w])
                else:
                    surfaces.add(int(w))
            tags = sorted(surfaces)
            c = field.add("Constant")
            field.setNumbers(c, "SurfacesList", tags)
            field.setNumber(c, "VIn", local)
            field.setNumber(c, "VOut", 1e22)
            field.setNumber(c, "IncludeBoundary", 1)
            self._fields.append(c)
            mins.append(c)
            if distance > 0 and size > local:
                d = field.add("Distance")
                field.setNumbers(d, "SurfacesList", tags)
                field.setNumber(d, "Sampling", 100)
                t = field.add("Threshold")
                field.setNumber(t, "InField", d)
                field.setNumber(t, "SizeMin", local)
                field.setNumber(t, "SizeMax", size)
                field.setNumber(t, "DistMin", 0.0)
                field.setNumber(t, "DistMax", distance)
                self._fields += [d, t]
                mins.append(t)
        m = field.add("Min")
        field.setNumbers(m, "FieldsList", mins)
        field.setAsBackgroundMesh(m)
        self._fields.append(m)

    def surface_mesh(self, orient: bool = True) -> SurfaceMesh:
        """The triangle surface mesh of the model as a :class:`SurfaceMesh`, for boundary element
        solvers. Raises :class:`MeshTopologyError` unless the mesh is watertight (every edge shared
        by exactly two triangles) and orientable. With ``orient`` the triangles are made
        consistent across surfaces and their normals point out of the solid (into the void for a
        cavity); the gmsh mesh is reversed accordingly, so :meth:`write` exports the same
        orientation. Call after ``mesh(dim=2)``; the triangle groups follow the face labels."""
        labels = self.face_tags() if self._faces else {}
        return build_surface_mesh(self.gmsh, labels, orient=orient)

    def write(self, path: str) -> None:
        """Write mesh/geometry; format from extension (.msh, .vtk, .stl, .step, .brep ...)."""
        self.gmsh.write(path)


def _interior_point(face: TrimmedSurface) -> tuple[float, float, float]:
    """A point inside the trimmed face (first grid point of the parameter domain it contains)."""
    s = face.surface
    (u0, u1), (v0, v1) = s.domain_u, s.domain_v
    for n in (8, 32, 128):
        for i in range(1, n):
            for j in range(1, n):
                u, v = u0 + (u1 - u0) * (i + 0.37) / n, v0 + (v1 - v0) * (j + 0.61) / n
                if u < u1 and v < v1 and face.contains(u, v):
                    return tuple(s(u, v))
    return tuple(s(0.5 * (u0 + u1), 0.5 * (v0 + v1)))


def _closest_surface(gmsh, tags: list[int], p) -> int:
    """The gmsh surface (of ``tags``) closest to the point ``p``."""
    return min(tags, key=lambda t: sum((a - b) ** 2 for a, b in
                                       zip(gmsh.model.getClosestPoint(2, t, p)[0], p)))


def add_shell(shell: Shell, solid: bool = True, physical_groups: bool = False,
              volume_name: str = "volume") -> list[int]:
    """Add a shell with its B-rep topology: faces share their edges. For a closed shell and
    ``solid=True`` a volume is created. Returns the gmsh surface tag of every face (in the order
    of ``shell.faces``, without repetitions) and the volume tag last, if one was created.

    The shell is handed over through Nurbforge's native STEP writer (exact NURBS, shared edges)
    and OCC's STEP reader, so OCC builds the topology itself – including the seam edges of closed
    surfaces (cylinders, cones, spheres, tori), which a face-by-face transfer with duplicate
    removal or sewing cannot reproduce once such a face is trimmed (e.g. a bore or a slot from a
    Boolean).

    ``physical_groups=True`` creates one physical surface group per face label
    (``TrimmedSurface.label``, e.g. set with ``Shell.with_face_labels`` before a Boolean;
    unlabelled faces go to ``"unlabeled"``) and a physical volume ``volume_name``, so boundary
    conditions can be assigned by name in the solver. Faces are matched to the gmsh surfaces by
    an interior point."""
    gmsh = _gmsh()
    occ = gmsh.model.occ
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "shell.step")
        write_step(path, [shell], unit="mm")
        unit = gmsh.option.getString("Geometry.OCCTargetUnit")
        gmsh.option.setString("Geometry.OCCTargetUnit", "MM")  # no unit conversion
        try:
            entities = occ.importShapes(path)
        finally:
            gmsh.option.setString("Geometry.OCCTargetUnit", unit)
    occ.synchronize()
    volumes = [t for d, t in entities if d == 3]
    if volumes:
        surfaces = sorted({t for _, t in gmsh.model.getBoundary([(3, v) for v in volumes],
                                                                 combined=False, oriented=False)})
    else:
        surfaces = [t for d, t in entities if d == 2]
    if volumes and not solid:
        occ.remove([(3, v) for v in volumes])  # keeps the surfaces
        occ.synchronize()
        volumes = []
    face_tags = [_closest_surface(gmsh, surfaces, _interior_point(f)) for f in shell.faces]
    if physical_groups:
        groups: dict[str, list[int]] = {}
        for f, tag in zip(shell.faces, face_tags):
            groups.setdefault(f.label or "unlabeled", []).append(tag)
        for name, group in groups.items():
            gmsh.model.addPhysicalGroup(2, sorted(set(group)), name=name)
        if volumes:
            gmsh.model.addPhysicalGroup(3, volumes, name=volume_name)
    return list(dict.fromkeys(face_tags)) + volumes


def export_step(objects, path: str, solid: bool = True) -> None:
    """Write surfaces, trimmed surfaces and shells to a STEP (or BREP/IGES, by extension) file
    via gmsh's OpenCASCADE kernel. Closed shells become solids if ``solid`` is set."""
    with GmshSession() as g:
        for o in objects:
            if isinstance(o, Shell):
                add_shell(o, solid=solid)
            else:
                g.add_surface(o)
        g.gmsh.model.occ.synchronize()
        g.write(path)


def remesh(mesh: SurfaceMesh, size: float, angle: float = 40.0) -> SurfaceMesh:
    """Remesh a closed triangle mesh with well-shaped triangles of edge length ``size``.

    gmsh splits the mesh into patches along edges where the normals of neighbouring triangles
    differ by more than ``angle`` (degrees), parametrises the patches and meshes them anew; the
    new nodes lie on the old triangles, so edges sharper than ``angle`` are kept and the shape
    changes by the chord error of the new triangles only. Returns a watertight,
    outward-oriented :class:`SurfaceMesh` with linear triangles (surface tags: the patches);
    raises :class:`MeshTopologyError` if the result is not closed. Face groups are not kept.
    Works inside an open :class:`GmshSession` as well, in a model of its own.

    The mesh must not contain nearly degenerate triangles (``ValueError`` if the worst
    :meth:`SurfaceMesh.quality` is below 0.01): their normals are noise, and gmsh's
    classification does not terminate on them. ``shrink_wrap`` returns meshes that qualify
    unless called with ``min_quality=0``."""
    if not size > 0.0:
        raise ValueError("remesh: size must be positive")
    if mesh.quality().min() < 0.01:
        raise ValueError("remesh: the mesh contains nearly degenerate triangles "
                         "(quality below 0.01)")
    gmsh = _gmsh()
    options = ("Mesh.MeshSizeMin", "Mesh.MeshSizeMax", "Mesh.MeshSizeFromCurvature",
               "Mesh.MeshSizeFromPoints", "Mesh.MeshSizeExtendFromBoundary",
               "Mesh.ElementOrder", "Mesh.RecombineAll")
    nested = bool(gmsh.isInitialized())
    if nested:
        previous = gmsh.model.getCurrent()
        saved = {name: gmsh.option.getNumber(name) for name in options}
    else:
        gmsh.initialize()
        gmsh.option.setNumber("General.Terminal", 0)
    try:
        gmsh.model.add("nurbforge-remesh")
        corners = mesh.triangles[:, :3]
        used = np.unique(corners)
        tags = np.arange(1, len(used) + 1)
        index = np.zeros(len(mesh.nodes), dtype=np.int64)
        index[used] = tags
        surface = gmsh.model.addDiscreteEntity(2)
        gmsh.model.mesh.addNodes(2, surface, tags.tolist(), mesh.nodes[used].ravel().tolist())
        gmsh.model.mesh.addElementsByType(surface, 2, [], index[corners].ravel().tolist())
        gmsh.model.mesh.classifySurfaces(math.radians(angle), True, True, math.pi)
        gmsh.model.mesh.createGeometry()
        for name, value in zip(options, (size, size, 0, 0, 0, 1, 0)):
            gmsh.option.setNumber(name, value)
        gmsh.model.mesh.generate(2)
        return build_surface_mesh(gmsh)
    finally:
        if nested:
            gmsh.model.remove()
            for name, value in saved.items():
                gmsh.option.setNumber(name, value)
            if previous:
                gmsh.model.setCurrent(previous)
        else:
            gmsh.finalize()


def bem_mesh(shell: Shell, size: float | None = None, order: int = 1,
             curvature: int | None = None, refine: dict | None = None,
             path: str | None = None) -> SurfaceMesh:
    """One-shot helper for boundary element solvers: mesh the closed ``shell`` with triangles and
    return a watertight, outward-oriented :class:`SurfaceMesh` (see
    :meth:`GmshSession.surface_mesh`). ``refine`` maps face labels or surface tags to local
    sizes; ``path`` additionally writes the mesh (.msh, .stl, .vtk, ...) with physical groups
    per face label."""
    with GmshSession() as g:
        g.add_shell(shell, solid=False, physical_groups=path is not None)
        for where, local in (refine or {}).items():
            g.refine(where, local)
        g.mesh(dim=2, size=size, order=order, curvature=curvature)
        mesh = g.surface_mesh()
        if path is not None:
            g.write(path)
    return mesh


def mesh_surfaces(surfaces: Iterable[NurbsSurface | TrimmedSurface], path: str,
                  size: float | None = None,
                  order: int = 1) -> None:
    """One-shot helper: add surfaces, mesh in 2D and write to ``path``."""
    with GmshSession() as g:
        for s in surfaces:
            g.add_surface(s)
        g.mesh(dim=2, size=size, order=order)
        g.write(path)
