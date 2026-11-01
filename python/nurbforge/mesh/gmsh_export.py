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
    """Context manager around ``gmsh.initialize()`` / ``gmsh.finalize()``."""

    def __init__(self, model_name: str = "nurbforge", verbose: bool = False):
        self.model_name = model_name
        self.verbose = verbose
        self.gmsh = _gmsh()

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

    def mesh(self, dim: int = 2, size: float | None = None, order: int = 1,
             heal: bool = True) -> None:
        """Synchronise the OCC model and generate a mesh."""
        g = self.gmsh
        if heal:
            # merges duplicate control-point vertices and fixes seams/poles
            g.model.occ.removeAllDuplicates()
        g.model.occ.synchronize()
        if size is not None:
            g.option.setNumber("Mesh.MeshSizeMax", size)
        g.model.mesh.generate(dim)
        if order > 1:
            g.model.mesh.setOrder(order)

    def write(self, path: str) -> None:
        """Write mesh/geometry; format from extension (.msh, .vtk, .stl, .step, .brep ...)."""
        self.gmsh.write(path)


def add_shell(shell: Shell, solid: bool = True) -> list[int]:
    """Add all faces of a shell; faces are merged along shared edges. For a closed shell and
    ``solid=True`` a volume is built from the surface loop. Returns the face tags (and the volume
    tag last, if one was created)."""
    gmsh = _gmsh()
    occ = gmsh.model.occ
    tags = [add_trimmed_surface(f) for f in shell.faces]
    occ.removeAllDuplicates()
    occ.synchronize()
    tags = [t for _, t in gmsh.model.getEntities(2)]
    if solid and shell.is_closed:
        loop = occ.addSurfaceLoop(tags, sewing=True)
        tags.append(occ.addVolume([loop]))
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
