# STEP export (native writer)

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
- One product, no assembly structure, colours, layers or names per face; no reader.
- gmsh's `importShapes` drops free faces when a file also contains solids (also for files written
  by OCC itself), so mixed files only show the solids there.
- Faces are written per shell; edges shared between different shells are written once per shell.
