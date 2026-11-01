# Nurbforge

[![CI](https://github.com/PhotonPower/Nurbforge/actions/workflows/ci.yml/badge.svg)](https://github.com/PhotonPower/Nurbforge/actions/workflows/ci.yml)
[![Wheels](https://github.com/PhotonPower/Nurbforge/actions/workflows/wheels.yml/badge.svg)](https://github.com/PhotonPower/Nurbforge/actions/workflows/wheels.yml)
[![Docs](https://github.com/PhotonPower/Nurbforge/actions/workflows/docs.yml/badge.svg)](https://photonpower.github.io/Nurbforge/)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

**Headless, NURBS-based CAD kernel for scientific software** – written in C++20 with a Python API.

Nurbforge builds exact 2D/3D curves, surfaces and surface systems for downstream numerical
tools such as optical design, FEM, BEM and FMM solvers. There is no GUI: geometry is created
from code (C++ or Python), from a local LLM (Ollama), and meshed with [gmsh](https://gmsh.info).

> **Status: early development (0.1.0).** The NURBS core, primitives, revolution/extrusion/ruled
> surfaces, affine transforms, arrays, gmsh transfer and the Ollama tool interface work and are
> tested. Intersections, trimming and Boolean operations are on the [roadmap](docs/roadmap.md).

## Features

| Area | Available now | Planned |
|---|---|---|
| NURBS core | Basis functions + derivatives, curve/surface evaluation, rational derivatives, knot insertion and refinement, Bézier decomposition, degree elevation and reduction, knot removal, clamping, splitting/joining, reversal, iso-curves | Interpolation, approximation |
| Primitives | Point, line, polyline, exact arc and circle, bilinear patch, rectangle | Ellipse |
| Optics | Exact conic surfaces, even/odd aspheres, Zernike (Fringe/Noll), XY polynomials, arbitrary sag functions – fitted to a tolerance | Trimmed circular apertures |
| Fitting | Curve/surface interpolation (incl. end tangents), least-squares approximation, fit to a tolerance | Constrained fits, smoothing |
| Construction | Extrusion, revolution (full/partial), ruled surface (auto degree/knot compatibility), loft/skinning, sweep (one/two rails, rotation-minimising frame, twist), Coons patches, Gordon surfaces, extrapolation (natural / linear), offset curves (planar) and surfaces with verified tolerance | Fillets |
| Transforms | Translate, rotate, scale, mirror, composition | – |
| Arrays | Linear, grid, polar | Path arrays |
| Analysis | Point, tangent, normal (incl. poles), partial derivatives of any order, global closest point / point inversion, curve–curve, curve–surface and surface–surface intersection, planar sections, curvature/torsion, Gaussian/mean/principal curvatures, arc length, area, sampling | B-rep topology, Booleans |
| Topology | Trimmed surfaces (parameter-space loops, exact circles), split along surface intersections, area by Green's theorem, exact gmsh export | B-rep, Boolean union/intersection/difference |
| Meshing | Exact NURBS transfer to gmsh (OCC), 2D meshing, high order, `.msh/.vtk/.stl/.step` | Physical groups, size fields, BEM-ready closed shells |
| LLM | Ollama tool-calling agent with whitelisted geometry tools, CLI | Constraint-aware planning, verification loop |

## Quick start

```bash
git clone https://github.com/PhotonPower/Nurbforge.git
cd Nurbforge
python -m pip install -e ".[dev]"      # builds the C++ core + Python bindings
pytest
```

```python
import math
import nurbforge as nf
from nurbforge.mesh import mesh_surfaces

# exact sphere: revolve a semicircle about z
meridian = nf.make_arc((0, 0, 0), (0, 0, -1), (1, 0, 0), 10.0, 0.0, math.pi)
sphere = nf.revolve(meridian, (0, 0, 0), (0, 0, 1), 2 * math.pi)

print(sphere(0.25, 0.5), sphere.normal(0.25, 0.5))
mesh_surfaces([sphere], "sphere.msh", size=1.0)
```

Natural language via a local model (requires a running `ollama serve`):

```bash
nurbforge-llm --model qwen2.5:14b --out result.msh
> Create a cylinder of radius 10 and height 30 along z
```

C++ only:

```bash
cmake --preset release && cmake --build --preset release && ctest --preset release
```

## Platforms

| Platform | Toolchain | Guide |
|---|---|---|
| AlmaLinux 10 | GCC 14, CMake ≥ 3.24, Python ≥ 3.10 | [docs/build/almalinux10.md](docs/build/almalinux10.md) |
| Windows 11 | Visual Studio 2022 or 2026 (MSVC), CMake ≥ 3.24 (≥ 4.2 for VS 2026), Python ≥ 3.10 | [docs/build/windows11.md](docs/build/windows11.md) |

## Binary wheels

Wheels for Linux (`manylinux_2_28` x86_64 – installs on AlmaLinux/RHEL 8, 9 and 10) and Windows
(`win_amd64`), CPython 3.10–3.13, are built by the [Wheels workflow](.github/workflows/wheels.yml)
for every `v*` tag and can be downloaded from its run artifacts:

```bash
pip install nurbforge-0.1.0-cp312-cp312-win_amd64.whl   # no compiler needed
```

## Documentation

The documentation site – these pages plus the Python (autodoc) and C++ (Doxygen/Breathe) API
references – is published at **<https://photonpower.github.io/Nurbforge/>**.

- [Architecture](docs/architecture.md) – layers, modules, data flow
- [NURBS fundamentals & conventions](docs/math/nurbs.md)
- [Python API overview](docs/api/python.md) · [C++ API overview](docs/api/cpp.md)
- [Optical surfaces](docs/optics.md)
- [Meshing with gmsh](docs/meshing.md)
- [LLM interface (Ollama)](docs/llm.md)
- [Roadmap](docs/roadmap.md) · [Architecture decisions](docs/adr/)
- [Contributing](CONTRIBUTING.md) · [Changelog](CHANGELOG.md)

## Repository layout

```
include/nurbforge/   public C++ headers (core, nurbs, geometry, modeling)
src/                 C++ implementation
python/bindings/     pybind11 module  -> nurbforge._core
python/nurbforge/    Python package (mesh/ = gmsh, llm/ = Ollama)
tests/cpp|python/    Catch2 and pytest suites
examples/            C++ and Python examples
docs/                documentation (Markdown), ADRs, Sphinx conf.py, Doxyfile
```

## License

MIT – see [LICENSE](LICENSE). Algorithms follow L. Piegl & W. Tiller, *The NURBS Book*, 2nd ed., Springer 1997.
