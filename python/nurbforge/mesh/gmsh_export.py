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
from collections.abc import Iterable

from .._core import NurbsCurve, NurbsSurface, Shell, TrimmedSurface, unique_knots


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
    """Add all faces of a shell; faces are merged along shared edges. For a closed shell and
    ``solid=True`` a volume is built from the surface loop. Returns the face tags (and the volume
    tag last, if one was created).

    ``physical_groups=True`` creates one physical surface group per face label
    (``TrimmedSurface.label``, e.g. set with ``Shell.with_face_labels`` before a Boolean;
    unlabelled faces go to ``"unlabeled"``) and a physical volume ``volume_name``, so boundary
    conditions can be assigned by name in the solver. Faces are matched to the gmsh surfaces
    (renumbered by the duplicate removal) by an interior point."""
    gmsh = _gmsh()
    occ = gmsh.model.occ
    for f in shell.faces:
        add_trimmed_surface(f)
    occ.removeAllDuplicates()
    occ.synchronize()
    tags = [t for _, t in gmsh.model.getEntities(2)]
    if physical_groups:
        groups: dict[str, list[int]] = {}
        for f in shell.faces:
            tag = _closest_surface(gmsh, tags, _interior_point(f))
            groups.setdefault(f.label or "unlabeled", []).append(tag)
    if solid and shell.is_closed:
        loop = occ.addSurfaceLoop(tags, sewing=True)
        tags.append(occ.addVolume([loop]))
        occ.synchronize()
    if physical_groups:
        for name, group in groups.items():
            gmsh.model.addPhysicalGroup(2, sorted(set(group)), name=name)
        if solid and shell.is_closed:
            gmsh.model.addPhysicalGroup(3, [tags[-1]], name=volume_name)
    return tags


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


def mesh_surfaces(surfaces: Iterable[NurbsSurface | TrimmedSurface], path: str,
                  size: float | None = None,
                  order: int = 1) -> None:
    """One-shot helper: add surfaces, mesh in 2D and write to ``path``."""
    with GmshSession() as g:
        for s in surfaces:
            g.add_surface(s)
        g.mesh(dim=2, size=size, order=order)
        g.write(path)
