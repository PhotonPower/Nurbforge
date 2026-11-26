# Nurbforge

[![CI](https://github.com/PhotonPower/Nurbforge/actions/workflows/ci.yml/badge.svg)](https://github.com/PhotonPower/Nurbforge/actions/workflows/ci.yml)
[![Wheels](https://github.com/PhotonPower/Nurbforge/actions/workflows/wheels.yml/badge.svg)](https://github.com/PhotonPower/Nurbforge/actions/workflows/wheels.yml)
[![Docs](https://github.com/PhotonPower/Nurbforge/actions/workflows/docs.yml/badge.svg)](https://photonpower.github.io/Nurbforge/)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

**Headless, NURBS-based CAD kernel for scientific software** – written in C++20 with a Python API.

Nurbforge builds exact 2D/3D curves, surfaces and surface systems for downstream numerical
tools such as optical design, FEM, BEM and FMM solvers. There is no GUI: geometry is created
from code (C++ or Python), from a local LLM (Ollama), and meshed with [gmsh](https://gmsh.info).

> **Status: 0.16.0 – early, API not yet stable.** NURBS core and toolbox (refinement, fitting, loft,
> sweep, Coons/Gordon, offsets, fillets, chamfers, curves on surfaces), freeform optical surfaces, intersections, trimming,
> B-rep shells, solid primitives and Boolean operations, solid operations for simulation (split,
> cap, thicken, hollow, pipes, defeaturing, edge fillets and chamfers), planar region Booleans, mass properties, shell
> sections and hatching (slicing for 2PP), gmsh transfer and BEM-ready surface meshes, STEP
> export and import, parallel use from Python (GIL-free operations, picklable objects) and the
> Ollama tool interface work and are tested on Linux, Windows and macOS. Next steps are on the [roadmap](docs/roadmap.md).

## Features

| Area | Available now | Planned |
|---|---|---|
| NURBS core | Basis functions + derivatives, curve/surface evaluation, rational derivatives, knot insertion and refinement, Bézier decomposition, degree elevation and reduction, knot removal, clamping, splitting/joining, reversal, iso-curves | – |
| Primitives | Point, line, polyline, exact arc and circle (also through 3 points), regular polygon, rectangle outline, helix / conical helix / spiral (C2, verified tolerance), bilinear patch, rectangle; closed solids: box, sphere, cylinder, cone/frustum, torus | Ellipse |
| Optics | Exact conic surfaces, even/odd aspheres, Zernike (Fringe/Noll), XY polynomials, arbitrary sag functions – fitted to a tolerance, circular apertures via `TrimmedSurface` | – |
| Fitting | Curve/surface interpolation (incl. end tangents), least-squares approximation, fit to a tolerance, best-fit line and plane (orthogonal regression) | Constrained fits, smoothing |
| Construction | Extrusion, revolution (full/partial), ruled surface (auto degree/knot compatibility), loft/skinning, sweep (one/two rails, rotation-minimising frame, twist), Coons patches, Gordon surfaces, extrapolation (natural / linear), offset curves (planar) and surfaces with verified tolerance, fillets (exact arcs between planar curves, all corners of a planar curve, rolling-ball fillets between surfaces with contact p-curves), chamfers (between curves, all corners of a curve), tween curves, tangent arcs and lines, exact projection onto planes, curves on surfaces (interpolation, geodesic offset), sweeps with several profiles, rail revolve, tapered extrusion, ribbons, cones to a point, three-sided patches, exact Coons patches from arcs, rebuild with reported error, surface blend/match (G0/G1/G2), two-rail sweeps with several profiles and non-uniform scaling, variable surface offset, chamfers between surfaces, developability check and unrolling of developable surfaces | Flattening of doubly curved surfaces |
| Transforms | Translate, rotate, scale, mirror, composition – of curves, surfaces and shells | – |
| Arrays | Linear, grid, polar | Path arrays |
| Analysis | Point, tangent, normal (incl. poles), partial derivatives of any order, global closest point / point inversion, curve–curve, curve–surface and surface–surface intersection, planar sections, sections of B-rep shells (oriented loops, critical heights, extrusion bands), hatching of section loops (exact hatch lines, inward contours that split, merge and vanish), curvature/torsion, Gaussian/mean/principal curvatures, arc length, area, sampling, curve analysis (osculating circle, curvature graph, direction and sense of rotation, G0/G1/G2 continuity between curves and at knots, duplicate curves), surface analysis (curvature and draft angle maps with extrema, direction, G0/G1/G2 continuity between surfaces and across shell edges, duplicate surfaces), silhouette curves (parallel and perspective views, cusps), hidden-line removal, 2D drawings with visible and hidden lines (Make2D), outlines | – |
| Topology | Trimmed surfaces (parameter-space loops, exact circles), split along surface intersections, area by Green's theorem, exact gmsh export; B-rep shells (sewing, watertightness, orientation, volume, centroid, inertia tensor), STEP export (native AP214 writer, or via gmsh/OCC) and import (native reader: B-spline and analytic faces), Boolean union/intersection/difference (incl. coincident faces and tangential contacts), point-in-solid; fillets and chamfers of solid edges (closed loops, convex and concave, corners of two or three edges, tangent chains, variable radius); shells with several components (voids, lumps) | Dedicated `Body` type |
| Meshing | Exact NURBS transfer to gmsh (OCC), 2D meshing, high order, `.msh/.vtk/.stl/.step`, physical groups from face labels (kept through Booleans, transforms and STEP), local mesh sizes by face label (with growth into the volume) and from curvature, BEM-ready surface meshes (watertight, outward-oriented, numpy arrays) | – |
| LLM | Ollama tool-calling agent with whitelisted geometry tools (incl. solids, mass properties, polygons, helices, rounded corners), CLI, verification rounds (`inspect` / `check`), constraint solving (`solve`) | – |

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

**Getting started:** the [tutorial](docs/tutorial.md) walks from curves to surfaces, trimmed
solids, Booleans, meshing and the LLM interface; every step is a runnable script in
[examples/python/](examples/python/).

C++ only:

```bash
cmake --preset release && cmake --build --preset release && ctest --preset release
```

## Platforms

| Platform | Toolchain | Guide |
|---|---|---|
| AlmaLinux 10 | GCC 14, CMake ≥ 3.24, Python ≥ 3.10 | [docs/build/almalinux10.md](docs/build/almalinux10.md) |
| Windows 11 | Visual Studio 2022 or 2026 (MSVC), CMake ≥ 3.24 (≥ 4.2 for VS 2026), Python ≥ 3.10 | [docs/build/windows11.md](docs/build/windows11.md) |
| macOS ≥ 14 (Apple Silicon, Intel) | Xcode Command Line Tools (Apple Clang), CMake ≥ 3.24, Python ≥ 3.10 | [docs/build/macos.md](docs/build/macos.md) |

## Binary wheels

Wheels for Linux (`manylinux_2_28` x86_64 – installs on AlmaLinux/RHEL 8, 9 and 10), Windows
(`win_amd64`) and macOS 14+ (`arm64`, `x86_64`), CPython 3.10–3.14, are built by the [Wheels workflow](.github/workflows/wheels.yml)
for every `v*` tag, published on [PyPI](https://pypi.org/project/nurbforge/) and attached, with the sdist, to the
[GitHub release](https://github.com/PhotonPower/Nurbforge/releases):

```bash
pip install nurbforge            # no compiler needed; nurbforge[mesh] adds gmsh
```

## Documentation

The documentation site – these pages plus the Python (autodoc) and C++ (Doxygen/Breathe) API
references – is published at **<https://photonpower.github.io/Nurbforge/>**.

- **User manual (PDF)** – the complete guide with a full command reference in the appendix;
  attached to every [GitHub release](https://github.com/PhotonPower/Nurbforge/releases) as
  `nurbforge-manual.pdf`, sources in [docs/manual](docs/manual/) (`python docs/manual/build.py`)
- **Reference for coding agents** – the whole API in one compact Markdown file for LLMs and
  agents, shipped with the wheel: `python -m nurbforge` prints it, `nurbforge.reference()`
  returns it ([python/nurbforge/REFERENCE.md](python/nurbforge/REFERENCE.md), also attached to
  every release; [llms.txt](llms.txt))
- [Tutorial](docs/tutorial.md) – the main workflow, step by step with runnable examples
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
