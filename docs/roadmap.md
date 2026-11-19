# Roadmap

## 0.1 – Foundation ✅
- NURBS core (basis, evaluation, rational derivatives, knot insertion)
- Primitives: line, polyline, arc, circle, bilinear patch, rectangle
- Extrusion, revolution, ruled surface; affine transforms; linear/grid/polar arrays
- Python bindings, gmsh transfer, Ollama agent; CI on AlmaLinux 10, Windows 11 and macOS

## 0.2 – Complete NURBS toolbox ✅
- Degree elevation ✅, degree reduction ✅, knot removal ✅
- Knot refinement ✅, Bézier decomposition ✅
- Curve/surface splitting and joining ✅
- Interpolation and least-squares approximation (curves and surfaces) ✅
- Full surface derivatives ✅, limit normals at poles ✅
- Closest point / point inversion ✅
- Arc length, curvature, torsion, area ✅
- Extrapolation (curve/surface extension: natural C∞ and linear G1/C1) ✅

## 0.3 – Surface construction
- Loft / skinning ✅
- Sweep (one and two rails) ✅
- Coons and Gordon patches ✅
- Offset curves/surfaces ✅
- Fillets ✅ (rolling ball): curve/curve arcs (exact), surface/surface rolling-ball fillets with verified tolerance; all corners of a planar curve (`fillet_corners`) ✅
- Curve constructors ✅ (#16): circle/arc through 3 points, regular polygon, rectangle outline, helix/conical helix/spiral with verified tolerance
- Optical surfaces: conics (exact), aspheres (even/odd), Zernike, XY polynomials, freeform sag functions → NURBS fit with verified tolerance ✅

## 0.4 – Intersections & trimming
- Curve–curve, curve–surface intersection ✅
- Surface–surface intersection (subdivision + marching) ✅
- Trimmed surfaces (parameter-space loops) ✅, split along intersections ✅, exact gmsh export ✅
- Surface systems with shared edges ✅ (sewn shells)
- Planar sections / cuts ✅

## 0.5 – Topology & Booleans
- B-rep shells (vertex/edge/coedge/face/shell) ✅, watertightness and orientation checks ✅, volume ✅; bodies with several shells ✅ (voids, lumps)
- Solid primitives as closed shells (box, sphere, cylinder, cone, torus) ✅, mass properties (centroid, inertia tensor, principal axes) ✅, affine transforms of shells ✅ (#18)
- Boolean union / intersection / difference on trimmed NURBS shells ✅ (incl. coincident faces and tangential contacts)
- STEP / IGES export via OCC ✅
- Native STEP writer (AP214, no dependencies) ✅
- Native STEP reader ✅ (B-spline and analytic geometry → exact NURBS, trim loops rebuilt in parameter space, units, voids)

## 0.6 – Slicing (for FourierGrating's CAD slicer) ✅
- Planar sections of shells as oriented closed loops (`section`), exact for extrusions / revolutions ✅
- Critical heights and extrusion bands along a direction ✅
- STEP: micrometre unit, solid names (`read_step_named`) ✅
- Thread safety of const queries documented and checked under ThreadSanitizer ✅

## 0.7 – Writing paths (for Voxelspaetzle's 2PP process model) ✅
- Hatching of section loops: exact hatch lines (holes, islands, tangential contacts) and inward contours ✅ (#44)
- Contours that split, merge and vanish (offset beyond the local half-width) ✅

## 0.8 – macOS and face labels ✅
- macOS (Apple Silicon and Intel, macOS ≥ 14): CI and binary wheels ✅
- Face labels through Booleans, transforms and STEP; gmsh physical groups ✅ (#18)

## 0.9 – Meshing for solvers ✅
- Local mesh sizes by face label or surface tag, growth with distance, curvature-based sizes ✅
- BEM-ready closed shells: watertight surface meshes with consistent outward orientation ✅

## 0.10 – Solids for simulation (#18) ✅
- Split a solid by a plane (symmetry), cut faces labelled ✅
- Cap planar openings of a shell (exact trim curves, nested openings as holes) ✅
- Thicken an open, smooth shell into a solid (offset + ruled side walls, symmetric option) ✅
- Hollow a smooth solid (inner void, or open at labelled faces) ✅
- Pipe along a path with constant, tapered or function radius, end faces labelled inlet/outlet ✅
- Defeaturing: remove small through and blind holes ✅, merge faces on the same surface ✅
- LLM tools for all solid operations ✅

## 0.11 – Parallel parameter studies ✅
- Expensive operations release the GIL in Python (Booleans, solid operations, intersections, sections, offsets, sweeps, fits, mass properties) ✅
- All geometry types picklable (processes, deepcopy) ✅
- Concurrent constructive operations checked under ThreadSanitizer ✅

## 0.12 – 2D profiles (#16) ✅
- Booleans of planar regions (exact curves, holes, touching regions), exact region area ✅
- Blend curves between curve ends (G1 cubic, G2 quintic) ✅
- Projection of curves onto surfaces (normal or directional, 3D curve + p-curve) ✅

## 0.13 – Curve tools (#16) ✅
- Chamfers between two curves and at all corners of a curve ✅
- Best-fit line and plane through points (orthogonal regression) ✅
- Exact parallel projection of curves and surfaces onto a plane ✅
- Tween curves; arcs from start, tangent and end; tangent lines from a point and common tangents of two curves ✅
- Curves on surfaces: interpolation through points in the parameter domain, geodesic offset within a surface ✅

## 0.14 – Edge chamfers and fillets on solids (#18) ✅
- Constant-radius fillets and chamfers of edges, edge chains and closed edge loops of a shell, convex (material removed) and concave (material added); adjacent faces cut back along the exact contact curves, end faces rerouted along the exact end curves ✅

## 0.15 – Vertex blends, variable radius, surface tools (#17, #18) ✅
- Edge chains meeting at a vertex: tangent continuation, mitre (mutual trimming), spherical / planar corner patch for three edges ✅
- Variable radius / chamfer distance along edge chains (smooth step between given points, exact arc sections) ✅
- Surface tools (#17): sweep with several profiles, rail revolve, tapered extrusion, ribbon, extrusion to a point, three-sided patch, exact Coons patches for rational boundaries, rebuild of curves and surfaces ✅
- Next from #17: surface blend / match (G1/G2), fitting to scattered data, variable offset, unrolling

## 0.16 – Surface blend / match, fitting to scattered data (#17)
- Blend surface between two surface sides and matching a surface to a target side, G0/G1/G2 ✅
- Surface fitting to scattered points (optionally with boundary curves and smoothing)

## Infrastructure
- ✅ Wheels via cibuildwheel (manylinux_2_28 for AlmaLinux 10 compatibility, win_amd64, macOS arm64/x86_64)
- ✅ Sphinx site with Breathe (Doxygen) API docs on GitHub Pages
- Benchmarks ✅ (`benchmarks/`, smoke-tested in CI); parallel seed search in the surface/surface intersection ✅ (2.3× on 16 threads); OpenMP on Boolean face pairs measured, not adopted (see benchmarks page)
- LLM: geometry verification step ✅ (`inspect` / `check` tools, verification rounds); constraint reasoning ✅ (`solve`: unknowns from measured targets)
