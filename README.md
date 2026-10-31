# Nurbforge

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
| Primitives | Point, line, polyline, exact arc and circle, bilinear patch, rectangle | Ellipse, conics, aspheres, Zernike/XY-polynomial optical surfaces |
| Construction | Extrusion, revolution (full/partial), ruled surface (auto degree/knot compatibility) | Loft/skinning, sweep, Coons/Gordon patches, offset, fillets |
| Transforms | Translate, rotate, scale, mirror, composition | – |
| Arrays | Linear, grid, polar | Path arrays |
| Analysis | Point, tangent, normal, sampling | Closest point, curvature, length/area, intersections |
| Topology | – | Trimming, B-rep, Boolean union/intersection/difference, extrapolation |
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

## Documentation

- [Architecture](docs/architecture.md) – layers, modules, data flow
- [NURBS fundamentals & conventions](docs/math/nurbs.md)
- [Python API overview](docs/api/python.md)
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
docs/                documentation, ADRs, Doxygen config
```

## License

MIT – see [LICENSE](LICENSE). Algorithms follow L. Piegl & W. Tiller, *The NURBS Book*, 2nd ed., Springer 1997.
