# Changelog

All notable changes are documented here. Format: [Keep a Changelog](https://keepachangelog.com),
versioning: [SemVer](https://semver.org).

## [Unreleased]
### Added
- Named tolerances `nf::tol::kLength`, `kParametric`, `kAngular` and the per-operation
  `nf::Tolerance` struct (`nurbforge/core/tolerance.hpp`, ADR 0005).
- Build guides for AlmaLinux 10 and Windows 11 (`docs/build/`).
- Knot refinement `NurbsCurve::refine_knots` (A5.4) and `NurbsSurface::refine_knots_u/v` (A5.5).
- Bézier decomposition `decompose_bezier()` for curves and surfaces; segments keep their original
  parameter intervals, unclamped knot vectors are supported.
- Degree elevation `NurbsCurve::elevate_degree` (A5.9) and `NurbsSurface::elevate_degree_u/v`.
- `NurbsCurve::clamped()` / `is_clamped()`: equivalent curve with clamped knot vector.
- Knot removal within a Cartesian tolerance (A5.8): `NurbsCurve::try_remove_knot` /
  `remove_knot`, `NurbsSurface::remove_knot_u/v`.
- Splitting and joining: `NurbsCurve::split` / `sub_curve`, `NurbsSurface::split_u/v` /
  `sub_surface` / `clamped`, free functions `join`, `join_u`, `join_v` (`nurbs/join.hpp`).
- `NurbsSurface::try_remove_knot_u/v`.
- Global closest point / point inversion on curves and surfaces (`analysis/closest_point.hpp`):
  branch and bound over Bézier pieces with safe lower bounds, safeguarded (projected) Newton;
  local variant with start value; Python `closest_point`, `closest_point_local`, `Tolerance`;
  LLM tool `closest_point`.
- Optical surfaces (`optics/surfaces.hpp`, docs/optics.md): exact conic surfaces of revolution,
  even/odd aspheres, Zernike (Fringe/Noll) and XY polynomial sags, generic `sag_surface`
  graph fit with adaptive knots and verified tolerance; `approximate_curve` with given knots
  and parameters (`fitting`); LLM tool `asphere`.
- Coons patches and Gordon surfaces (`modeling/gordon.hpp`) as Boolean sums with automatic
  boundary orientation, network-consistency check and non-rational conversion of rational
  curves; LLM tool `coons`.
- Sweep (`modeling/sweep.hpp`): exact translational sweep, sweep along a trajectory with
  rotation-minimising frame (double reflection) and optional twist, two-rail sweep with uniform
  scaling; adaptive piecewise skinning to a tolerance; LLM tool `sweep`.
- Loft / skinning through section curves (`modeling/loft.hpp`), incl. rational sections;
  `make_compatible` for any number of curves; interpolation with explicit parameters and
  averaged data parameters (`fitting`); LLM tool `loft`.
- Extrapolation (`modeling/extend.hpp`): natural (C∞, conics stay exact) and linear (C1/G1)
  extension of curves at either end and of surfaces across any side by an arc length, keeping
  the original parameters; LLM tool `extend`.
- Fitting (`fitting/fit.hpp`): curve interpolation (optionally with end derivatives),
  least-squares approximation, approximation to a tolerance, surface interpolation and
  approximation of point grids; uniform/chord-length/centripetal parameters; LLM tool
  `interpolate_points`.
- Curve curvature, torsion, arc length and its inverse; surface Gaussian/mean/principal
  curvatures with principal directions (also at poles) and area (`analysis/properties.hpp`,
  adaptive Gauss–Kronrod); LLM tool `measure`.
- `NurbsSurface::derivatives(u, v, order)`: all partial derivatives of any order (A3.6 + A4.4).
- Degree reduction with rigorous error verification: `NurbsCurve::reduce_degree`,
  `NurbsSurface::reduce_degree_u/v`.
### Changed
- All tolerance literals in `src/` replaced by the named constants; angular comparisons now
  uniformly use 1e-9.
### Fixed
- `NurbsSurface::normal` returns the limit normal at poles and collapsed edges instead of throwing
  (e.g. at the poles of revolved spheres – relevant for BEM meshes).
- `ruled` / `make_compatible` now handle curves of different degree (previously
  `NotImplementedError`) and unclamped curves.
- `NurbsCurve::is_rational` compares weights relatively, so uniformly scaled weights
  (e.g. all 1e-12) are no longer misclassified.
- `pip install .` on Windows with Visual Studio 2026 and an older system CMake: CMake >= 4.2 is
  now required there and fetched automatically.

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
