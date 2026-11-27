# Meshing with gmsh

`nurbforge.mesh` passes NURBS data **exactly** into gmsh's OpenCASCADE kernel
(`occ.addBSpline`, `occ.addBSplineSurface`) – degree, distinct knots, multiplicities and weights.
gmsh therefore meshes the true geometry. The test suite verifies that gmsh's own evaluation
(`gmsh.model.getValue`) agrees with Nurbforge to 1e-9.

```python
from nurbforge.mesh import GmshSession

with GmshSession(verbose=False) as g:
    g.add_surface(front)
    g.add_surface(back)
    g.mesh(dim=2, size=0.5, order=2)   # second-order elements, e.g. for BEM
    g.write("lens.msh")                # .msh, .vtk, .stl, .step, .brep, ...
```

Inside the session the full gmsh API (`g.gmsh`) is available for anything not wrapped here.

Trimmed surfaces are transferred exactly as well - the trim loops become B-spline wires in the
parameter space of the surface:

```python
aperture = nf.TrimmedSurface(zernike_surface, nf.circle_loop(0.0, 0.0, 10.0))  # u = x, v = y
with GmshSession() as g:
    g.add_surface(aperture)        # TrimmedSurface or NurbsSurface
    g.mesh(dim=2, size=0.5)
```

Mesh nodes on trimmed edges follow OCC's own 3D approximation of the trim curve
(deviation ~1e-4 in our tests).

Shells keep their B-rep topology (shared edges, seam edges of closed faces); closed shells become
solids:

```python
from nurbforge.mesh import export_step
shell = nf.Shell.sew([nf.TrimmedSurface(side), bottom, top]).oriented()
export_step([shell], "cylinder.step")      # .step / .brep / .iges by extension
```

Note: OCC's `getMass` is only accurate to ~0.5 % for such faces; Nurbforge's own `area()` /
`volume()` are accurate to the requested tolerance (the volume of a fine tetrahedral mesh of the
re-imported solid agrees).

Without gmsh, `nf.write_step(path, shells)` writes STEP files with the kernel's native writer
(shared edges kept, units selectable) – see [io.md](io.md).

## Physical groups from face labels

Faces carry a label (`TrimmedSurface.label`). The solid primitives label their faces (box
`bottom`/`top`/`front`/`back`/`left`/`right`, cylinder and cone `side`/`bottom`/`top`),
`Shell.with_face_labels` sets them, and Booleans, transforms and STEP keep them – every piece of
a split face inherits the label of its face. `add_shell(shell, physical_groups=True)` turns them
into named physical surface groups (unlabelled faces: `"unlabeled"`) plus a physical volume, so
boundary conditions are assigned by name after any chain of operations:

```python
import nurbforge as nf
from nurbforge.mesh import GmshSession, add_shell

pipe = nf.make_cylinder((0, 0, 0), (0, 0, 1), 1.0, 5.0).with_face_labels(["wall", "inlet", "outlet"])
part = nf.boolean(pipe, nf.make_box((-2, -0.2, 2), (4, 0.4, 1)), nf.BooleanOp.Difference)
with GmshSession() as g:
    g.add_shell(part, physical_groups=True)   # wall, inlet, outlet, bottom/top/... of the slot, volume
    g.mesh(dim=3, size=0.3)
    g.write("pipe.msh")
```

`g.face_tags()` returns the current gmsh surface tags of every label, e.g. for gmsh calls of your
own. (The module function `add_shell` works too, but only the session method remembers the
labels for `face_tags` and `refine`.)

## Local mesh sizes

`refine` sets a mesh size on faces, given by label (shells added with `g.add_shell`) or gmsh
surface tag. With `distance` the size grows linearly from the local size on the faces to the
global size at that distance – also into the volume, so the elements near an inlet or a
contact face stay small without refining the whole part. `curvature` makes elements smaller on
strongly curved faces (number of elements per 2π):

```python
pipe = nf.make_cylinder((0, 0, 0), (0, 0, 1), 1.0, 5.0).with_face_labels(["wall", "inlet", "outlet"])
with GmshSession() as g:
    g.add_shell(pipe, physical_groups=True)
    g.refine("inlet", 0.05, distance=1.0)        # 0.05 on the inlet, 0.3 from 1.0 away
    g.refine("outlet", 0.1)                      # only the outlet face itself
    g.mesh(dim=3, size=0.3, curvature=24)        # global size, 24 elements per 2π of curvature
    g.write("pipe.msh")
```

Several refinements combine by taking the smallest size. The faces get their size exactly
(gmsh `Constant` field); the growth uses a `Distance` + `Threshold` field, and the background
field is the minimum of all. Faces adjacent to a refined face also become finer near the shared
edges, because gmsh grades the mesh away from the fine boundary. Labels are resolved to gmsh
surfaces when meshing, so `refine` may be called before or after further `add_shell` calls and
duplicate removal (which renumbers surfaces).

## Surface meshes for boundary element solvers

BEM solvers need a closed triangle mesh whose normals point consistently out of the body.
`bem_mesh` meshes a closed shell and returns such a mesh as numpy arrays:

```python
from nurbforge.mesh import bem_mesh

part = nf.boolean(nf.make_box((-2, -2, -2), (4, 4, 4)), nf.make_sphere((0, 0, 0), 1.0),
                  nf.BooleanOp.Difference)           # box with a spherical cavity
m = bem_mesh(part, size=0.2, order=2, refine={"surface": 0.1}, path="part.msh")
m.nodes          # (N, 3) coordinates
m.triangles      # (M, 3) or (M, 6) for order=2 – 0-based, gmsh node order
m.groups         # triangle indices per face label, e.g. m.groups["top"]
m.normals(), m.areas(), m.centroids(), m.volume()   # volume() == part.volume() up to the mesh error
```

In a session, `g.surface_mesh()` does the same after `g.mesh(dim=2)`. It checks the mesh and
raises `MeshTopologyError` if an edge is used by one triangle (open shell, or surfaces added
separately with `add_surface` that do not share their edges) or by more than two
(non-manifold), or if the mesh contains quadrangles. It then orients the mesh:

- the triangles of neighbouring surfaces are made consistent (every shared edge is traversed in
  opposite directions);
- every closed component is turned so that its normals point **out of the solid** – away from
  the material on the outer boundary and into the void on a cavity (decided by its signed
  volume and how many other components enclose it).

The orientation is applied to the gmsh model as well (`gmsh.model.mesh.reverse`), so
`g.write("part.msh")` / `.stl` / `.vtk` exports the same, BEM-ready mesh. `orient=False`
returns gmsh's mesh unchanged (still checked for watertightness); `m.is_watertight()`,
`m.is_consistently_oriented()` and `m.boundary_edges()` check any `SurfaceMesh`.

The volume of a quadratic mesh is computed from its corner triangles, so it converges like a
linear mesh.

## Shrink wrap: a closed hull around anything

`shrink_wrap` builds one closed triangle mesh around shells and loose surfaces that overlap,
leave gaps or are open – a fallback where an exact Boolean or sewing is not possible, a repair
for imported geometry, and a way to get a simplified outer boundary. Unlike the rest of
Nurbforge it is **approximate**: the result is a mesh accurate to about the `resolution`, not
exact geometry. It does not need gmsh: the work is done by the C++ kernel (`nf::shrink_wrap` in
`mesh/shrink_wrap.hpp`, see the [C++ overview](api/cpp.md)).

```python
from nurbforge.mesh import shrink_wrap

hull = shrink_wrap([body_a, body_b, *loose_surfaces],
                   resolution=0.2,     # grid cell: about the edge length and the accuracy
                   offset=0.5,         # distance of the hull from the input (default 0)
                   fill_holes=True)
hull.is_watertight(), hull.volume()    # a SurfaceMesh, normals pointing outwards
hull.write("envelope.msh")             # or .stl
```

- `objects`: a `Shell`, `TrimmedSurface` or `NurbsSurface`, or a list of them.
- Gaps and openings narrower than about 1.5 x `resolution` (plus twice the `offset`) are
  closed; wider ones are not, and the hull then follows the inside of an open part as well.
  Parts further apart give separate closed components.
- With `offset=0` the hull lies on the input: on a sphere the nodes are within about
  0.1 x `resolution` of the surface, sharp edges are bevelled at the size of the grid, and
  concave corners are rounded with a radius of 0.75 x `resolution`.
- Parts thinner than about 1.5 x `resolution` vanish at `offset=0`. A loose surface has no
  thickness: give it an `offset` of at least 0.75 x `resolution`. If nothing is left,
  `ValueError` is raised.
- `fill_holes=True` (default) fills everything that cannot be reached from outside, so only
  the outer hull remains. With `False`, enclosed voids stay as inner components with normals
  pointing into the void – including the inside of a closed shell, since the input is treated
  as a set of surfaces; the result is then a wall of thickness 2 x `offset`.
- The grid covers the bounding box; `max_nodes` (default 20 million) limits its size, and
  `ValueError` tells the size it would have had.

Method: the input is sampled into points a quarter of the resolution apart; the distance to
them is computed on the grid in a narrow band; the nodes reachable from the border without
coming closer to the input than the hull distance are outside; the surface between them and
the rest is extracted with marching tetrahedra, which always gives a closed, consistently
oriented 2-manifold. That hull is watertight only at a distance of at least
0.75 x `resolution`, so for a smaller `offset` a second distance field, measured from the
hull, gives the surface that lies the difference further in (a morphological closing).

### Triangle quality and remeshing

The triangles cut out of the grid can be arbitrarily thin. `m.quality()` gives the shape
quality of every triangle of a `SurfaceMesh`: 4·√3·area / (sum of the squared sides), 1 for an
equilateral triangle and 0 for a degenerate one. `shrink_wrap` checks it: if the worst triangle
is below `min_quality` (default 0.3), the nodes are moved within the surface towards the
centres of their neighbours until it is reached (at most 50 smoothing steps). The connectivity
is kept, so the mesh stays closed, and since the nodes move in their tangent planes the shape
and the volume change by fractions of a per cent of the resolution only. Typically 5 steps
give a worst quality above 0.3 and a mean of 0.85. `min_quality=0` returns the mesh as cut out
of the grid.

For triangles of a different size, or of still better shape (mean quality above 0.95),
`remesh` builds a new mesh through gmsh:

```python
from nurbforge.mesh import remesh

coarse = remesh(hull, size=0.5)        # angle=40: edges sharper than 40 degrees are kept
```

gmsh splits the mesh into patches along its sharp edges, parametrises them and meshes them
anew; the new nodes lie on the old triangles. The result is again a watertight,
outward-oriented `SurfaceMesh`. `remesh` refuses meshes with nearly degenerate triangles
(quality below 0.01, `ValueError`), on which gmsh's classification does not terminate.

## Implementation notes

- Control-point indexing differs: Nurbforge is row-major with v fastest; gmsh expects u fastest.
  `add_surface` reorders accordingly.
- Control points are created as temporary gmsh points and removed after the spline is built;
  otherwise they would remain as free vertices and appear as stray mesh nodes.
- `mesh(heal=True)` calls `occ.removeAllDuplicates()` to merge coincident edges of adjacent patches.
- `add_shell` hands the shell over through Nurbforge's native STEP writer and OCC's STEP reader
  (a temporary file, exact NURBS), so OCC builds the topology itself. A face-by-face transfer
  cannot reproduce the seam of a *trimmed* closed face (cylinder with a bore or slot, sphere,
  cone or torus cut by a Boolean): the outer trim loop runs along the seam twice, duplicate
  removal then deletes such faces and sewing leaves broken volumes.

## Installation

`pip install gmsh` (or `pip install nurbforge[mesh]`). On AlmaLinux 10 the pip wheel needs a few
X/GL runtime libraries even when used headless – see [build/almalinux10.md](build/almalinux10.md).
