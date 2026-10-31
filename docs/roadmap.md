# Roadmap

## 0.1 – Foundation ✅
- NURBS core (basis, evaluation, rational derivatives, knot insertion)
- Primitives: line, polyline, arc, circle, bilinear patch, rectangle
- Extrusion, revolution, ruled surface; affine transforms; linear/grid/polar arrays
- Python bindings, gmsh transfer, Ollama agent; CI on AlmaLinux 10 and Windows 11

## 0.2 – Complete NURBS toolbox
- Degree elevation ✅, degree reduction ✅, knot removal ✅
- Knot refinement ✅, Bézier decomposition ✅
- Curve/surface splitting and joining ✅
- Interpolation and least-squares approximation (curves and surfaces)
- Full surface derivatives ✅, limit normals at poles ✅
- Closest point / point inversion ✅
- Arc length, curvature, torsion, area ✅
- Extrapolation (curve/surface extension with G1/G2 continuity)

## 0.3 – Surface construction
- Loft / skinning, sweep (one and two rails), Coons and Gordon patches
- Offset curves/surfaces, fillets (rolling ball, approximate)
- Optical surfaces: conics/aspheres (even/odd), Zernike, XY polynomials, freeforms → NURBS fit with error bound

## 0.4 – Intersections & trimming
- Curve–curve, curve–surface, surface–surface intersection (subdivision + marching)
- Trimmed surfaces (parameter-space loops), surface systems with shared edges
- Planar sections / cuts

## 0.5 – Topology & Booleans
- B-rep (vertex/edge/loop/face/shell/body), watertightness checks
- Boolean union / intersection / difference on trimmed NURBS shells
- STEP / IGES export (via OCC or native writer)

## Infrastructure
- Wheels via cibuildwheel (manylinux_2_28 for AlmaLinux 10 compatibility, win_amd64)
- Sphinx site with Breathe (Doxygen) API docs on GitHub Pages
- Benchmarks; optional OpenMP parallel sampling
- LLM: geometry verification step (agent checks results), constraint reasoning
