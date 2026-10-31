# Changelog

All notable changes are documented here. Format: [Keep a Changelog](https://keepachangelog.com),
versioning: [SemVer](https://semver.org).

## [Unreleased]

## [0.1.0] – 2026-10-31
### Added
- C++20 NURBS kernel: knot vectors, basis functions and derivatives, rational curves and surfaces,
  rational derivatives, knot insertion, reversal, iso-curves, sampling.
- Primitives: line, polyline, exact arc/circle, bilinear patch, rectangle.
- Construction: extrusion, revolution, ruled surface with automatic knot compatibility.
- Affine transforms (translate, rotate, scale, mirror) and linear/grid/polar patterns.
- pybind11 Python package `nurbforge` with NumPy interop.
- Exact gmsh transfer and meshing (`nurbforge.mesh`).
- Ollama tool-calling agent and CLI `nurbforge-llm` (`nurbforge.llm`).
- Catch2 and pytest suites, CI for AlmaLinux 10 and Windows 11, documentation and ADRs.
