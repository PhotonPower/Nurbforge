# Architecture

## Layers

```
 ┌──────────────────────────────────────────────────────────────┐
 │  Python package  nurbforge                                   │
 │   ├─ llm/   Ollama client, whitelisted tools, agent, CLI     │
 │   ├─ mesh/  exact NURBS transfer to gmsh (OCC), meshing      │
 │   └─ patterns.py, __init__ (re-exports of _core)             │
 ├──────────────────────────────────────────────────────────────┤
 │  nurbforge._core  (pybind11, python/bindings/module.cpp)     │
 ├──────────────────────────────────────────────────────────────┤
 │  C++ kernel  libnurbforge_core  (namespace nf)               │
 │   topology/  trimmed surfaces, split, shells (sewing)        │
 │   optics/    conics, aspheres, Zernike, XY polynomials       │
 │   modeling/  construct, loft, sweep, Coons/Gordon, extend,   │
 │              offset, fillet                                  │
 │   fitting/   interpolation, least-squares approximation      │
 │   analysis/  closest point, curvature, arc length, area      │
 │   geometry/  transform, primitives                           │
 │   nurbs/     knot_vector, basis, curve, surface, join,       │
 │              refinement / elevation / reduction / splitting  │
 │   core/      Vec3/Vec4, tolerance, errors                    │
 └──────────────────────────────────────────────────────────────┘
```

Dependencies point strictly downwards. The C++ kernel has **no third-party runtime
dependencies**; gmsh and Ollama are integrated on the Python side only, so the kernel stays
embeddable into solver codes.

## Design principles

1. **Exactness first.** Conics, spheres, cylinders and tori are represented exactly as rational
   NURBS. Data handed to gmsh is the original NURBS (knots, multiplicities, weights), never a
   tessellation.
2. **Value semantics, immutability.** Geometry objects are immutable values; operations return
   new objects (`curve.transformed(t)`, `curve.insert_knot(u)`). This makes them thread-safe and
   easy to expose to Python.
3. **Validated invariants.** Constructors validate knot vectors, weights and net sizes and throw
   `nf::GeometryError` (Python: `nurbforge.GeometryError`, a `ValueError`).
4. **Headless & scriptable.** Every operation is reachable from Python; the LLM layer only
   composes existing API calls.
5. **Portable.** C++20 standard library only; tested with GCC (AlmaLinux 10) and MSVC (Windows 11).

## Data model

| Type | Storage |
|---|---|
| `Vec3` / `Vec4` | Cartesian / homogeneous points |
| `NurbsCurve` | degree p, knot vector U (size n+p+2), control points P_i, weights w_i |
| `NurbsSurface` | degrees (p,q), knots U,V, net `num_u × num_v`, row-major `i*num_v + j` (i along u) |
| `Transform` | affine 3×4 matrix, row-major |

`TrimmedSurface`: surface + closed parameter-space loops (outer CCW, holes CW).
`Shell`: trimmed faces sewn along shared edges (vertices, edges, coedges).
Planned: `Body` (several shells, solids),
tolerance context per model.

## Error handling

- `GeometryError` – invalid input or degenerate geometry.
- `NotImplementedError` – recognised request for a planned feature (raised by the Python layer
  and reserved in the kernel for roadmap items that are not implemented yet).

## Build

CMake (≥ 3.24) builds the static core, the tests (Catch2, fetched if absent), examples and the
pybind11 module. `pip install .` uses scikit-build-core with the same CMake files.
