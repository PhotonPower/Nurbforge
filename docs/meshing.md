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

Inside the session the full gmsh API (`g.gmsh`) is available, e.g. for physical groups or size fields.

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

Shells are added face by face with duplicate removal; closed shells become solids:

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
    add_shell(part, physical_groups=True)     # wall, inlet, outlet, bottom/top/... of the slot, volume
    g.mesh(dim=3, size=0.3)
    g.write("pipe.msh")
```

## Implementation notes

- Control-point indexing differs: Nurbforge is row-major with v fastest; gmsh expects u fastest.
  `add_surface` reorders accordingly.
- Control points are created as temporary gmsh points and removed after the spline is built;
  otherwise they would remain as free vertices and appear as stray mesh nodes.
- `mesh(heal=True)` calls `occ.removeAllDuplicates()` to merge coincident edges of adjacent patches.

## Installation

`pip install gmsh` (or `pip install nurbforge[mesh]`). On AlmaLinux 10 the pip wheel needs a few
X/GL runtime libraries even when used headless – see [build/almalinux10.md](build/almalinux10.md).
