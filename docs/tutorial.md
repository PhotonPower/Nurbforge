# Tutorial

This tutorial walks through the main Nurbforge workflow in Python: curves → surfaces → fitted
and freeform surfaces → intersections and trimming → closed solids and Booleans → offsets and
fillets → meshing → the LLM interface. Every snippet is taken from a runnable script in
[`examples/python/`](https://github.com/PhotonPower/Nurbforge/tree/main/examples/python); run
the scripts to see the printed numbers (each one finishes in a few seconds). The test
`tests/python/test_examples.py` runs all of them and checks that the snippets below still match
the scripts.

```bash
python -m pip install -e ".[dev]"    # builds the C++ core and the Python bindings
python examples/python/01_primitives.py
```

Background: the [NURBS conventions](math/nurbs.md) page defines the parametrisations and
explains each algorithm, the [Python API overview](api/python.md) lists all functions, and the
[Python API reference](api/python-reference.md) has the full signatures.

## Conventions in a nutshell

- Points and vectors are `nf.Vec3`; anywhere a point is expected you can pass a tuple or list.
  NumPy rows are *not* converted implicitly – use `tuple(row)` or `nf.Vec3(*row)`. Bulk inputs
  (point lists and grids for fitting) take NumPy arrays of shape `(n, 3)` / `(nu, nv, 3)`.
- Curves are `nf.NurbsCurve`, evaluated with `c(u)`; surfaces are `nf.NurbsSurface`, evaluated
  with `s(u, v)`. `c.domain`, `s.domain_u`, `s.domain_v` give the parameter ranges.
- Invalid geometry raises `nf.GeometryError` (a `ValueError`) with a precise message.
- Approximating operations (fits, sweeps, offsets, optical sags, intersections) take a
  `tolerance` and verify it at dense check points; exact constructions (arcs, circles,
  revolutions, extrusions, conics) are exact up to rounding. See
  [ADR 0005](adr/0005-tolerance-model.md) for the tolerance model.

## 1. Curves and primitives

Lines, arcs and circles are exact (circles and arcs are rational quadratic NURBS):

```python
# from examples/python/01_primitives.py
line = nf.make_line((0, 0, 0), (10, 0, 0))
circle = nf.make_circle(center=(0, 0, 0), normal=(0, 0, 1), radius=5.0)
arc = nf.make_arc((0, 0, 0), (1, 0, 0), (0, 1, 0), 5.0, 0.0, math.radians(120))
```

`make_arc(center, x_axis, y_axis, radius, start, end)` runs from angle `start` to `end`
measured from `x_axis` towards `y_axis`. Copies are made with transforms and patterns:

```python
# from examples/python/03_patterns.py
hole = nf.make_circle((40, 0, 0), (0, 0, 1), 4.0)
holes = nf.apply_pattern(hole, nf.polar_pattern((0, 0, 0), (0, 0, 1), 8, 2 * math.pi))
```

## 2. Surfaces from curves

Revolution and extrusion are exact. A lens is three revolved surfaces
([02_lens_and_mesh.py](https://github.com/PhotonPower/Nurbforge/blob/main/examples/python/02_lens_and_mesh.py)):

```python
# from examples/python/02_lens_and_mesh.py
def spherical_cap(R, z_vertex, sign):
    """Cap with vertex at z_vertex; sign=+1 curves towards -z, -1 towards +z."""
    alpha = math.asin(h / R)
    center = (0, 0, z_vertex - sign * R)
    meridian = nf.make_arc(center, (0, 0, sign), (1, 0, 0), R, 0.0, alpha)
    return nf.revolve(meridian, (0, 0, 0), (0, 0, 1), 2 * math.pi)
```

Free-form surfaces through several curves: **loft** (skinning) through sections, **sweep** of a
profile along a path (rotation-minimising frame by default), **Coons** patches from four
boundary curves and **Gordon** surfaces through a curve network
([06_loft_sweep_gordon.py](https://github.com/PhotonPower/Nurbforge/blob/main/examples/python/06_loft_sweep_gordon.py)):

```python
# from examples/python/06_loft_sweep_gordon.py
sections = [nf.make_circle((0, 0, z), (0, 0, 1), r) for z, r in radii]
vase = nf.loft(sections, degree_v=3)  # u along the sections, v across them
```

```python
# from examples/python/06_loft_sweep_gordon.py
start, tangent = helix(0.0), helix.tangent(0.0)
profile = nf.make_circle(tuple(start), tuple(tangent), 0.15)  # profile at the path start
tube = nf.sweep(profile, helix, nf.SweepOptions(frame=nf.SweepFrame.RotationMinimizing,
                                                tolerance=1e-5))
```

The profile must be given at the start of the path. The tube stays within about 1e-6 of radius
0.15 around the helix. `coons(c0, c1, d0, d1)` needs the curves to form a closed loop with
`c0`/`c1` at v = 0 / 1 and `d0`/`d1` at u = 0 / 1 (curve directions are fixed automatically):

```python
# from examples/python/06_loft_sweep_gordon.py
patch = nf.coons(c0, c1, d0, d1)
```

Math: [Loft](math/nurbs.md#loft-skinning), [Sweep](math/nurbs.md#sweep),
[Coons and Gordon](math/nurbs.md#coons-patches-and-gordon-surfaces).

## 3. Fitting measured data

Interpolation passes through every point; least squares uses fewer control points;
`approximate_curve_tol` finds the smallest control net meeting a tolerance
([05_fitting.py](https://github.com/PhotonPower/Nurbforge/blob/main/examples/python/05_fitting.py)):

```python
# from examples/python/05_fitting.py
interp = nf.interpolate_curve(pts, degree=3)
lsq = nf.approximate_curve(pts, degree=3, n_ctrl=8)
fit = nf.approximate_curve_tol(pts, degree=3, tolerance=1e-4)
```

For 40 points on a spiral the script reports a deviation of ~1e-16 for the 40-point
interpolant, 8e-4 with 8 control points and 7e-5 with the 12 control points chosen for a
tolerance of 1e-4. Surfaces are fitted to point grids of shape `(nu, nv, 3)`:

```python
# from examples/python/05_fitting.py
grid = np.stack([x, y, np.sin(x) * np.cos(y)], axis=2)  # shape (nu, nv, 3)

surf = nf.interpolate_surface(grid, degree_u=3, degree_v=3)
coarse = nf.approximate_surface(grid, degree_u=3, degree_v=3, ctrl_u=6, ctrl_v=5)
```

`nf.data_parameters(points)` returns the parameters assigned to the data points (centripetal
by default), useful to measure the deviation at the data. Math: [Fitting](math/nurbs.md#fitting).

## 4. Optical and freeform surfaces

Aspheres, Zernike and XY-polynomial sags are fitted to a tolerance; pure conics are exact
([07_freeform_optics.py](https://github.com/PhotonPower/Nurbforge/blob/main/examples/python/07_freeform_optics.py)):

```python
# from examples/python/07_freeform_optics.py
asph = nf.Asphere(curvature=1 / 40.0, conic=-0.8, even=[1e-6, -2e-9])
r_max = 15.0
a_surf = nf.asphere_surface(asph, r_max, tolerance=1e-8)
```

```python
# from examples/python/07_freeform_optics.py
zern = nf.ZernikeSag(base=nf.Asphere(curvature=1 / 80.0), norm_radius=10.0,
                     indexing=nf.ZernikeIndexing.Fringe,
                     coeffs=[0, 0, 0, 0.002, 0.001, -0.0015, 0, 0, 0.0005])  # Z1..Z9
z_surf = nf.zernike_surface(zern, tolerance=1e-8)  # graph over [-10, 10]^2
```

The Zernike and XY surfaces are graphs with u = x, v = y, so `z_surf(x, y).z` is directly
comparable with `zern.sag(x, y)` (the script finds errors of a few 1e-9 for a tolerance of
1e-8). A circular clear aperture is a trimmed surface:

```python
# from examples/python/07_freeform_optics.py
aperture = nf.TrimmedSurface(z_surf, nf.circle_loop(0.0, 0.0, 9.0))
```

`nf.surface_curvature(s, u, v)` returns Gaussian, mean and principal curvatures, e.g. 1/R at
the vertex of a sphere. Details: [Optical surfaces](optics.md).

## 5. Intersections, trimming and sections

`intersect` works for curve/curve, curve/surface, curve/plane and surface/surface. Surface
intersections return 3D curves plus p-curves on both surfaces
([08_intersect_trim.py](https://github.com/PhotonPower/Nurbforge/blob/main/examples/python/08_intersect_trim.py)):

```python
# from examples/python/08_intersect_trim.py
res = nf.intersect(ball, cylinder)  # SurfaceIntersection: .curves and isolated .points
```

For a sphere of radius 2 and a cylinder of radius 1 the result is two circles at z = ±√3, each
returned in pieces that break at the seams of the closed surfaces (total length 2π to 1e-8).

Trim loops live in the (u, v) parameter space of the surface. On a 2 × 2 square (u = x/2,
v = y/2) a hole of radius 0.25 in (u, v) is a circle of radius 0.5:

```python
# from examples/python/08_intersect_trim.py
square = nf.make_rectangle((0, 0, 0), (2, 0, 0), (0, 2, 0))  # u = x / 2, v = y / 2
holed = nf.TrimmedSurface(square).with_hole(nf.circle_loop(0.5, 0.5, 0.25))  # r = 0.5 in xy
```

`holed.area()` gives 4 − π/4 exactly (Green's theorem on the exact trim circle). Splitting a
surface along its intersection with another yields trimmed regions, and planar sections slice a
surface with parallel planes:

```python
# from examples/python/08_intersect_trim.py
parts = nf.split(plane, ball)  # trimmed regions: the disk inside and the frame outside
```

```python
# from examples/python/08_intersect_trim.py
slices = nf.planar_sections(ball, nf.Plane((0, 0, -1.6), (0, 0, 1)), spacing=0.8, count=5)
```

Math: [Intersections](math/nurbs.md#intersections),
[Trimmed surfaces](math/nurbs.md#trimmed-surfaces).

## 6. Closed solids and Boolean operations

`Shell.sew` stitches trimmed faces along shared edges; `oriented()` makes the normals point
outwards; closed shells report their volume
([09_solids_booleans.py](https://github.com/PhotonPower/Nurbforge/blob/main/examples/python/09_solids_booleans.py)):

```python
# from examples/python/09_solids_booleans.py
def ball(center, r):
    """Sphere from a revolved semicircle: one face, closed by its seam edge and two poles."""
    meridian = nf.make_arc(center, (0, 0, -1), (1, 0, 0), r, 0.0, math.pi)
    return nf.Shell.sew([nf.TrimmedSurface(nf.revolve(meridian, center, (0, 0, 1), 2 * math.pi))])
```

A cylinder solid is the extruded circle plus two square lids trimmed to the inscribed circle:

```python
# from examples/python/09_solids_booleans.py
def cylinder(r, z0, z1):
    """Cylinder along z: extruded circle plus two lids trimmed to circles."""
    side = nf.extrude(nf.make_circle((0, 0, z0), (0, 0, 1), r), (0, 0, z1 - z0))
    lids = [nf.TrimmedSurface(nf.make_rectangle((-r, -r, z), (2 * r, 0, 0), (0, 2 * r, 0)),
                              nf.circle_loop(0.5, 0.5, 0.5)) for z in (z0, z1)]
    return nf.Shell.sew([nf.TrimmedSurface(side), *lids]).oriented()
```

Booleans take two closed shells and return a closed shell:

```python
# from examples/python/09_solids_booleans.py
r = 0.3
drilled = nf.boolean(cube, cylinder(r, -2.0, 2.0), nf.BooleanOp.Difference)
```

The script compares union, intersection and difference of two unit spheres at distance 1 with
the analytic lens volume π(4R + d)(2R − d)²/12 (relative errors ≤ 1e-7, limited by the fitted
intersection curves and the requested integration tolerance) and gets 8 − 2π·0.3² for the
drilled cube. `nf.point_in_shell(shell, p)` classifies points. Math:
[Shells](math/nurbs.md#shells-b-rep), [Boolean operations](math/nurbs.md#boolean-operations).

## 7. Offsets and fillets

Offsets keep the parameterisation and are verified against the requested tolerance; they raise
`GeometryError` where the distance reaches a radius of curvature
([10_offsets_fillets.py](https://github.com/PhotonPower/Nurbforge/blob/main/examples/python/10_offsets_fillets.py)):

```python
# from examples/python/10_offsets_fillets.py
d = -0.2  # d > 0 offsets to the left of the curve direction seen from plane_normal (inwards)
off = nf.offset_curve(loop, d, plane_normal=(0, 0, 1), tolerance=1e-7)
```

```python
# from examples/python/10_offsets_fillets.py
shell = nf.offset_surface(sphere, sign * 0.5, tolerance=1e-6)
```

The surface offset goes along `S_u × S_v`; check its direction with `surface.normal(u, v)` as
the script does. A curve fillet is an exact arc and comes with the trimmed input curves, ready to
be joined into a G1 chain; the rolling-ball fillet between two surfaces returns the fillet
surface, the spine and the contact curves with their p-curves on both surfaces:

```python
# from examples/python/10_offsets_fillets.py
f = nf.fillet_curves(a, b, 0.5)  # automatic side: the inner corner
```

```python
# from examples/python/10_offsets_fillets.py
chain = nf.join(nf.join(f.trimmed_a, f.arc, 1e-9), f.trimmed_b, 1e-9)  # G1 chain
```

```python
# from examples/python/10_offsets_fillets.py
sf = nf.fillet_surfaces(floor, wall, 0.5)  # ball on the +normal side of both: x > 0, z > 0
```

Math: [Offsets](math/nurbs.md#offsets), [Fillets](math/nurbs.md#fillets).

## 8. Meshing and export with gmsh

`nurbforge.mesh` hands the exact NURBS data (including trimmed surfaces and shells) to gmsh's
OpenCASCADE kernel. gmsh is optional; the examples catch its absence:

```python
# from examples/python/02_lens_and_mesh.py
try:
    mesh_surfaces([front, back, edge], "lens.msh", size=1.0)
    print("wrote lens.msh")
except (ImportError, OSError) as exc:  # gmsh missing or its runtime libraries unavailable
    print(f"meshing skipped: {exc}")
```

Closed shells become solids in STEP/BREP/IGES files:

```python
# from examples/python/09_solids_booleans.py
try:  # optional STEP export through gmsh's OpenCASCADE kernel
    from nurbforge.mesh import export_step

    export_step([drilled], "drilled_cube.step")
    print("wrote drilled_cube.step")
except (ImportError, OSError) as exc:  # gmsh missing or its runtime libraries unavailable
    print(f"STEP export skipped: {exc}")
```

For full control (element order, size fields, physical groups) use `GmshSession`, see
[Meshing with gmsh](meshing.md).

## 9. Geometry from natural language

A local model served by [Ollama](https://ollama.com) can call a whitelisted set of geometry
tools (it never executes code). This example needs a running `ollama serve`:

```python
# from examples/python/04_llm_ollama.py
agent = GeometryAgent(OllamaClient(model="qwen2.5:14b"))
print(agent.ask("Create a cylinder with radius 10 and height 30 along z, "
                "then a polar array of 6 copies of a circle of radius 2 at x=20 in the xy-plane."))
```

The created objects live in `agent.workspace` as ordinary Nurbforge objects. The same tools can
be called without a model, which is how the test suite exercises them:

```python
from nurbforge.llm import Workspace
from nurbforge.llm.tools import call_tool

ws = Workspace()
call_tool(ws, "make_circle", {"name": "c", "center": [0, 0, 0], "normal": [0, 0, 1], "radius": 1})
print(call_tool(ws, "closest_point", {"name": "c", "point": [3, 0, 0]}))  # ... distance=2 ...
```

Details and the tool list: [LLM interface](llm.md).

## Next steps

- C++: the same kernel without Python – see `examples/cpp/lens_surface.cpp`,
  `examples/cpp/booleans.cpp` and the [C++ API overview](api/cpp.md).
- [Architecture](architecture.md) for the module layout, [Roadmap](roadmap.md) for what is
  planned.
