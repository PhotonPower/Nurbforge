# STEP export and import (native writer and reader)

Nurbforge writes B-rep shells to STEP (ISO 10303-21) files with its own, dependency-free writer
in the C++ core (`include/nurbforge/io/step.hpp`). Unlike the gmsh route
(`nurbforge.mesh.export_step`, see [meshing.md](meshing.md)) it needs neither gmsh nor
OpenCASCADE, writes the kernel's B-spline data unchanged and keeps the shell's topology
(shared edges and vertices) as it is.

```python
import nurbforge as nf

drilled = nf.boolean(cube, cylinder, nf.BooleanOp.Difference)
nf.write_step("drilled_cube.step", drilled)               # a Shell or a list of Shells
nf.write_step("parts.step", [a, b], unit="m", name="parts", author="me")
text = nf.to_step(drilled)                                # the file contents as a string
```

```cpp
#include "nurbforge/io/step.hpp"
nf::write_step("drilled_cube.step", {drilled});
std::string text = nf::to_step({drilled}, {.unit = nf::StepLengthUnit::Metre});
```

## Options

| Option | Default | Meaning |
|---|---|---|
| `unit` | `"mm"` | Length unit declared in the file: `"mm"`, `"cm"`, `"m"`, `"inch"` (C++: `StepLengthUnit`). Coordinates are written unchanged; the unit says what one model unit means. Readers convert on import (OCC to mm: a 2 × 2 × 2 cube written with `unit="m"` imports with volume 8·10⁹). |
| `name` | `"nurbforge"` | `PRODUCT` name and `FILE_NAME` name |
| `author`, `organization` | `""` | `FILE_NAME` fields |
| `timestamp` | current UTC time | ISO 8601 time stamp in `FILE_NAME` (set it for reproducible files) |
| `uncertainty` | `1e-6` | `UNCERTAINTY_MEASURE_WITH_UNIT` (model units); OCC uses it as initial vertex / edge tolerance |

## File structure

- Header: `FILE_DESCRIPTION`, `FILE_NAME`, `FILE_SCHEMA(('AUTOMOTIVE_DESIGN { 1 0 10303 214 1 1 1 1 }'))`
  (AP214).
- Product structure: `APPLICATION_CONTEXT`, `APPLICATION_PROTOCOL_DEFINITION`, `PRODUCT`,
  `PRODUCT_DEFINITION_FORMATION`, `PRODUCT_DEFINITION`, `PRODUCT_DEFINITION_SHAPE`,
  `SHAPE_DEFINITION_REPRESENTATION`, one product for all shells.
- Representation context: `GEOMETRIC_REPRESENTATION_CONTEXT(3)` with the length unit
  (`SI_UNIT`, or `CONVERSION_BASED_UNIT('INCH', …)`), radian, steradian and the uncertainty.
- Closed shells: `MANIFOLD_SOLID_BREP(CLOSED_SHELL(…))` in an
  `ADVANCED_BREP_SHAPE_REPRESENTATION`. A shell that is not consistently oriented is oriented
  first (`Shell::oriented`). Every connected component is classified by its signed volume:
  components with outward normals are solids; components with inward normals (voids, e.g. a
  cube minus an inner sphere) become `ORIENTED_CLOSED_SHELL(…, .F.)` voids of the smallest
  solid containing them (`BREP_WITH_VOIDS`).
- Open shells: `SHELL_BASED_SURFACE_MODEL(OPEN_SHELL(…))` in a
  `MANIFOLD_SURFACE_SHAPE_REPRESENTATION`. If a call contains both kinds, solids and surface
  model go into one plain `SHAPE_REPRESENTATION`. Empty shells are skipped.
- Faces: `ADVANCED_FACE(bounds, surface, same_sense)`. The outer trim loop becomes the
  `FACE_OUTER_BOUND`, holes `FACE_BOUND`s, each an `EDGE_LOOP` of `ORIENTED_EDGE`s. Orientation
  is carried by the surface parametrisation (`TrimmedSurface::flipped` reparametrises), so
  `same_sense` is `.T.` except where a face must point against its surface normal (void shells,
  inside-out components): then `same_sense` and the bound orientations are `.F.`.
- Edges: every `ShellEdge` is one `EDGE_CURVE(VERTEX_POINT, VERTEX_POINT, curve, .T.)`, referenced
  by the `ORIENTED_EDGE`s of all its uses (orientation `.F.` where the loop runs against the
  curve). A seam of a closed surface appears twice in the same loop with opposite orientations.
  Pole points of revolved surfaces (degenerate boundary pieces) are not edges. Edges split at
  T-junctions while sewing are ordered along their loop by chaining their vertices.
- Geometry: `B_SPLINE_SURFACE_WITH_KNOTS` / `B_SPLINE_CURVE_WITH_KNOTS` with distinct knots and
  multiplicities (`knot_spec .UNSPECIFIED.`, `self_intersect .F.`, `closed` flags computed from
  the control net). Rational geometry uses the complex instances

  ```
  (BOUNDED_SURFACE() B_SPLINE_SURFACE(pu,pv,net,.UNSPECIFIED.,uc,vc,.F.)
   B_SPLINE_SURFACE_WITH_KNOTS(mu,mv,ku,kv,.UNSPECIFIED.) GEOMETRIC_REPRESENTATION_ITEM()
   RATIONAL_B_SPLINE_SURFACE(weights) REPRESENTATION_ITEM('') SURFACE())
  (BOUNDED_CURVE() B_SPLINE_CURVE(p,points,.UNSPECIFIED.,closed,.F.)
   B_SPLINE_CURVE_WITH_KNOTS(m,k,.UNSPECIFIED.) CURVE() GEOMETRIC_REPRESENTATION_ITEM()
   RATIONAL_B_SPLINE_CURVE(weights) REPRESENTATION_ITEM(''))
  ```

  (partial entity types in alphabetical order, each with only its own explicit attributes).
  Control points and weights are listed with the u index outermost, as stored.
- Numbers: shortest round-trip representation (`std::to_chars`), always with a decimal point
  (`1.`, `0.25`, `1.E-05`). Long instances are wrapped after commas.

## Validation

`tests/python/test_step.py` writes shells natively, re-imports them with OpenCASCADE (via gmsh,
if installed) and compares OCC's volume with `Shell.volume()`:

| Shape | `Shell.volume()` | OCC volume after import | Relative difference |
|---|---|---|---|
| Cube 2 × 2 × 2 | 8 | 8.000000000000 | 2·10⁻¹⁶ |
| Unit sphere (rational, seam) | 4.188790205 | 4.187243221 | 3.7·10⁻⁴ |
| Union of two unit spheres (distance 1) | 7.068583460 | 7.067161594 | 2.0·10⁻⁴ |
| Cube minus cylinder r = 0.3 (drilled) | 7.434513326 | 7.435302412 | 1.1·10⁻⁴ |
| Cube minus sphere r = 0.5 (void) | 7.476401224 | 7.476577524 | 2.4·10⁻⁵ |

The deviations for curved faces come from OCC's mass integration (planar faces are exact);
every import yields exactly one solid. `tests/cpp/test_step.cpp` checks the file syntax:
entity counts, every reference resolves, each edge is used twice with opposite orientations,
seams, complex-entity order and real-number syntax.

## Limitations

- Edge curves are the shell's 3D edge curves, i.e. the trim curves mapped to 3D and fitted
  (cubic, within 1/100 of the sewing tolerance), not exact conics; no p-curves
  (`SURFACE_CURVE` / `PCURVE`) are written – readers project the edges onto the faces (OCC does
  this without complaint).
- One product, no assembly structure, colours, layers or names per face.
- gmsh's `importShapes` drops free faces when a file also contains solids (also for files written
  by OCC itself), so mixed files only show the solids there.
- Faces are written per shell; edges shared between different shells are written once per shell.

## Reading STEP

The native reader (`src/io/step_read.cpp`) imports the B-rep solids and surface models of STEP
files written by other CAD systems (AP203 / AP214 / AP242 – the schema name is not checked) or
by `write_step`, without OpenCASCADE:

```python
shells = nf.read_step("part.step")                # list of Shells, converted to millimetres
shells = nf.read_step("part.step", unit="m")      # or "cm", "inch"
shells = nf.from_step(text, tolerance=1e-5)       # file contents; sewing tolerance (target units)
shells = nf.read_step("part.step", skip_unsupported=True)  # skip bad faces with a UserWarning
```

```cpp
#include "nurbforge/io/step.hpp"
std::vector<nf::Shell> shells = nf::read_step("part.step");
std::vector<std::string> warnings;
auto s = nf::from_step(text, {.unit = nf::StepLengthUnit::Metre, .skip_unsupported = true,
                              .warnings = &warnings});
```

Each `MANIFOLD_SOLID_BREP` becomes one shell, a `BREP_WITH_VOIDS` one shell holding the outer
boundary and the voids (void faces point into the void, as `Shell::oriented` expects, so
`volume()` is the solid volume), and every shell of a `SHELL_BASED_SURFACE_MODEL` one shell.
The faces are sewn with `Shell::sew(faces, tolerance)`. The LLM tool `import_step` stores the
shells as `name`, `name_1`, ….

### Supported entities

| Kind | Entities |
|---|---|
| Part 21 syntax | header and `DATA` sections, simple and complex instances `#n=(A(…) B(…) …)`, strings with `''`, enumerations, `$`, `*`, typed parameters (`LENGTH_MEASURE(1.)`), lists, comments, line breaks anywhere |
| Topology | `MANIFOLD_SOLID_BREP`, `BREP_WITH_VOIDS`, `SHELL_BASED_SURFACE_MODEL`, `CLOSED_SHELL`, `OPEN_SHELL`, `ORIENTED_CLOSED_SHELL`, `ORIENTED_OPEN_SHELL`, `ADVANCED_FACE` / `FACE_SURFACE` (`same_sense`), `FACE_OUTER_BOUND`, `FACE_BOUND` (orientation), `EDGE_LOOP`, `VERTEX_LOOP` (ignored), `POLY_LOOP`, `ORIENTED_EDGE`, `EDGE_CURVE` (`same_sense`), `VERTEX_POINT` |
| Surfaces | `B_SPLINE_SURFACE_WITH_KNOTS`, `BEZIER_SURFACE`, `UNIFORM_SURFACE`, `QUASI_UNIFORM_SURFACE` (also rational complex instances); `PLANE`, `CYLINDRICAL_SURFACE`, `CONICAL_SURFACE`, `SPHERICAL_SURFACE`, `TOROIDAL_SURFACE` |
| Curves | `B_SPLINE_CURVE_WITH_KNOTS`, `BEZIER_CURVE`, `UNIFORM_CURVE`, `QUASI_UNIFORM_CURVE` (also rational), `LINE`, `CIRCLE`, `ELLIPSE`, `HYPERBOLA`, `PARABOLA`, `POLYLINE`; `SURFACE_CURVE`, `SEAM_CURVE`, `INTERSECTION_CURVE` (their 3D curve), `TRIMMED_CURVE` (basis curve; the vertices bound the edge) |
| Placement | `CARTESIAN_POINT`, `DIRECTION`, `VECTOR`, `AXIS2_PLACEMENT_3D` |
| Units | `SI_UNIT` with any prefix, `CONVERSION_BASED_UNIT` (e.g. `'INCH'` via `LENGTH_MEASURE_WITH_UNIT`) for length and plane angle (cone semi-angles), from `GLOBAL_UNIT_ASSIGNED_CONTEXT`; files without units are read as millimetres |

Other entities met while translating a face (e.g. `OFFSET_SURFACE`, `SURFACE_OF_REVOLUTION`)
raise `GeometryError` naming the entity and instance (`unsupported entity OFFSET_SURFACE (#30)
(face #37)`); with `skip_unsupported` the face is left out (the shell stays open) and the
message is reported as a warning. `FACETED_BREP` is reported the same way. Product structure,
presentation (colours, layers), PMI and p-curves are ignored.

### Geometry translation

- B-spline surfaces and curves are taken as they are (clamped if the file uses unclamped knots).
- Analytic surfaces become exact rational NURBS patches whose normal `S_u × S_v` is the STEP
  surface normal: planes a rectangle around the face boundary (+5 %), cylinders and cones a
  revolved line, spheres a revolved half circle, tori a revolved circle. The angular range is
  chosen from the boundary: the largest angular gap that no boundary curve crosses and from which
  the boundary turns away (judged by the inward boundary normal) lies outside the face; without
  such a gap the full turn starts at a seam edge of the face (an edge used twice), so seams of
  the file and of the patch coincide. Cones stop at the apex.
- Edges are built from their vertices: lines and conic arcs exactly (circles / ellipses as
  rational arcs, hyperbolas and parabolas as rational quadratic Bézier segments), B-spline
  edges as the sub-curve between the vertex parameters (closed curves wrapped around their
  start); end points are set exactly to the vertices.

### Trim loops

Every edge use is mapped into the parameter space of its face surface:

- planar patches with affine parametrisation (the writer's rectangles, the patches of `PLANE`):
  exactly, by the inverse affine map applied to the edge curve;
- otherwise by projecting samples (local Newton from the previous sample, global closest point
  as fallback) and interpolating them with the edge parameter; the image of the interpolant is
  verified against the 3D edge between the samples and refined until it lies within
  `tolerance / 10` (plus the edge's own distance from the surface).

Seams of closed surfaces: samples on the seam take the side of their nearest regular sample;
curves lying on the seam (a seam edge appears twice in its loop) take the side given by the
orientation (the face lies to the left of the loop, seen against the face normal). Curves are
split where they cross the seam or pass through a pole. Points at a pole take the free
parameter of their neighbours, and consecutive curves that end at different parameters of the
same pole are connected along it. Loops that run once around a closed surface (a band on a
cylinder or torus without seam edge, a spherical cap) are joined along the seam or closed
through the pole. The loop with the largest parameter-space area becomes the outer loop.
Faces with `same_sense .F.` (and faces of `ORIENTED_CLOSED_SHELL(.F.)`) are flipped with
`TrimmedSurface::flipped`.

### Validation

`tests/cpp/test_step_read.cpp` reads back the writer's output (cube, sphere, cylinder with lids,
union of two spheres, cube with spherical void, open shell, several shells, metre and inch
files; face counts, closedness, orientation, area and volume within 10⁻⁶) and hand-written
files in the style of other CAD systems (cube from `PLANE` / `LINE`, cylinder from
`CYLINDRICAL_SURFACE` / `CIRCLE`, half cylinder, sphere with a `VERTEX_LOOP`, hemisphere).
`tests/python/test_step.py` reads files written by gmsh / OpenCASCADE (if gmsh is installed);
volumes after import (`Shell.volume(1e-7)`):

| OCC model | Volume read | Reference | Relative difference |
|---|---|---|---|
| Box 10 × 20 × 30 | 6000 | 6000 (exact) | 0 |
| Cylinder r = 5, h = 20 | 1570.796327 | 500 π | 10⁻¹⁶ |
| Cone frustum (seam, two circles) | 1455.604596 | 1390 π / 3 | 3·10⁻¹³ |
| Cone with apex | 261.799388 | 250 π / 3 | 5·10⁻¹² |
| Sphere r = 4 (`VERTEX_LOOP` only) | 268.0825731 | 256 π / 3 | 2·10⁻¹⁶ |
| Torus R = 20, r = 5 (seams in both directions) | 9869.604401 | 1000 π² | 2·10⁻¹³ |
| Box minus cylinder | 717.256661 | 1000 − 90 π | 10⁻¹⁵ |
| Box with all edges filleted (r = 2) | 907.7049927 | OCC 907.7049927 | 10⁻¹³ |
| Cone cut parallel to its axis (`HYPERBOLA` edges) | 65.6912391 | OCC 65.6912385 | 9·10⁻⁹ |
| Sphere cut through the poles | 56.5486678 | 18 π | 10⁻¹² |
| Box with spherical void (`BREP_WITH_VOIDS`) | 966.489678 | 1000 − 32 π / 3 | 10⁻¹⁵ |

### Limitations

- Assemblies: placement transformations (`MAPPED_ITEM`, `ITEM_DEFINED_TRANSFORMATION`,
  `CONTEXT_DEPENDENT_SHAPE_REPRESENTATION`) are not applied – all solids are returned in the
  coordinates of their own representation. One length unit per file (the first
  `GLOBAL_UNIT_ASSIGNED_CONTEXT`).
- Not supported: `SURFACE_OF_REVOLUTION`, `SURFACE_OF_LINEAR_EXTRUSION`, `OFFSET_SURFACE`,
  `RECTANGULAR_TRIMMED_SURFACE`, `CURVE_BOUNDED_SURFACE`, faceted B-reps, tori with minor radius
  ≥ major radius, wireframe and mesh geometry.
- P-curves in the file are not used (trim loops are always rebuilt from the 3D edges), so
  edges must lie on their faces within the tolerance; sloppy files may need a larger
  `tolerance`.
- On a closed B-spline surface, a face whose boundary crosses the seam twice without a seam edge
  (a patch straddling the seam) is rejected; the surface is not re-parametrised.
