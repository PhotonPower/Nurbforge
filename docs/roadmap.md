# Roadmap

## 0.1 – Foundation ✅
- NURBS core (basis, evaluation, rational derivatives, knot insertion)
- Primitives: line, polyline, arc, circle, bilinear patch, rectangle
- Extrusion, revolution, ruled surface; affine transforms; linear/grid/polar arrays
- Python bindings, gmsh transfer, Ollama agent; CI on AlmaLinux 10 and Windows 11

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
- Fillets (rolling ball, approximate)
- Optical surfaces: conics (exact), aspheres (even/odd), Zernike, XY polynomials, freeform sag functions → NURBS fit with verified tolerance ✅

## 0.4 – Intersections & trimming
- Curve–curve, curve–surface intersection ✅
- Surface–surface intersection (subdivision + marching) ✅
- Trimmed surfaces (parameter-space loops) ✅, split along intersections ✅, exact gmsh export ✅
- Surface systems with shared edges ✅ (sewn shells)
- Planar sections / cuts ✅

## 0.5 – Topology & Booleans
- B-rep shells (vertex/edge/coedge/face/shell) ✅, watertightness and orientation checks ✅, volume ✅; bodies with several shells
- Boolean union / intersection / difference on trimmed NURBS shells
- STEP / IGES export via OCC ✅ (native writer open)

## Infrastructure
- ✅ Wheels via cibuildwheel (manylinux_2_28 for AlmaLinux 10 compatibility, win_amd64)
- ✅ Sphinx site with Breathe (Doxygen) API docs on GitHub Pages
- Benchmarks; optional OpenMP parallel sampling
- LLM: geometry verification step (agent checks results), constraint reasoning
