# Changelog

All notable changes are documented here. Format: [Keep a Changelog](https://keepachangelog.com),
versioning: [SemVer](https://semver.org).

## [Unreleased]
### Added
- Native STEP writer (`io/step.hpp`, `docs/io.md`): `write_step(path, shells, options)` and
  `to_step(shells, options)` write ISO 10303-21 files with the AP214 schema without any
  dependency – closed shells as `MANIFOLD_SOLID_BREP` (`BREP_WITH_VOIDS` for inner void
  components), open shells as `SHELL_BASED_SURFACE_MODEL`; `ADVANCED_FACE`s on exact
  (rational) B-spline surfaces, one `EDGE_CURVE` per shell edge shared by its oriented uses
  (seams twice per loop), length unit mm / cm / m / inch, full-precision reals. Python
  `nf.write_step(path, shell_or_shells, unit=...)`, `nf.to_step`; LLM tool `export_step`.
  Round trips through OpenCASCADE (gmsh) give one solid each, volumes within 4·10⁻⁴ of
  `Shell.volume()` (exact for planar faces).
- Booleans with coincident faces: face pairs are classified by sampling (general / coincident /
  touching only), coincident faces are cut along each other's projected loops instead of being
  intersected, and common parts are classified "on" (same / opposite normal) – boxes sharing
  faces, stacked or identical solids, coaxial cylinders. Redundant (duplicate) cuts are dropped
  before splitting.
- Documentation: tutorial (`docs/tutorial.md`, in the Sphinx guide and linked from the README)
  walking from curves and surfaces through fitting, freeform optics, intersections and trimming,
  closed shells and Booleans, offsets and fillets to meshing and the LLM interface; new Python
  examples `05_fitting.py`, `06_loft_sweep_gordon.py`, `07_freeform_optics.py`,
  `08_intersect_trim.py`, `09_solids_booleans.py`, `10_offsets_fillets.py` (each prints its key
  numbers against analytic references) and the C++ example `examples/cpp/booleans.cpp`;
  `tests/python/test_examples.py` runs every example (also with gmsh hidden) and checks that the
  tutorial snippets are verbatim excerpts of the scripts. `02_lens_and_mesh.py` now runs without
  gmsh and prints sags and areas; the C++ examples build with the project warning flags.
- Fillets (`modeling/fillet.hpp`): `fillet_curves` (exact rational arc of radius r between two
  coplanar intersecting curves; centre from the offset curves refined by Newton, tangent
  parameters, trimmed curves forming a G1 chain; automatic or explicit sides) and
  `fillet_surfaces` (rolling-ball fillet: spine from the intersection of the offset surfaces,
  exact ball positions in the planes normal to the spine, exact circular sections skinned with
  verified tolerance; returns surface, spine, contact curves and their p-curves on both
  surfaces); Python `fillet_curves`, `fillet_surfaces`, `CurveFillet`, `SurfaceFillet`; LLM
  tool `fillet`.
- Offset curves and surfaces (`modeling/offset.hpp`): `offset_curve` (planar, cubic B-spline
  fitted per Bézier segment and joined, cusp/corner/planarity checks) and `offset_surface`
  (along the unit normal incl. limit normals at poles, quintic tensor-product fit with the knot
  lines of the input as break lines, curvature and crease checks); both keep the parameter
  correspondence and verify the tolerance at dense check points; Python `offset_curve`,
  `offset_surface`; LLM tool `offset`.
- Binary wheels via cibuildwheel (`.github/workflows/wheels.yml`, `[tool.cibuildwheel]` in
  `pyproject.toml`): manylinux_2_28 x86_64 and win_amd64, CPython 3.10–3.13, each wheel tested
  with the pytest suite; sdist; built on `v*` tags, manually and on packaging PRs. PyPI
  publishing is prepared (trusted publishing) but disabled.
- Sphinx documentation site (`docs/conf.py`): furo theme, MyST for the existing Markdown pages,
  autodoc Python reference, Breathe C++ reference from Doxygen XML (`build/docs/xml`); built by
  `.github/workflows/docs.yml` and deployed to GitHub Pages from `main`.
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
- Boolean operations on closed shells (`topology/boolean.hpp`): face/face intersection clipped to
  both trimmed faces, general face splitting by a planar arrangement in the parameter domain
  (`split(TrimmedSurface, cuts)`), point-in-solid classification by ray casting; LLM tool
  `boolean`. Sewing resolves T-junctions.
- B-rep shells (`topology/shell.hpp`): sewing of trimmed faces along shared edges, free and
  non-manifold edges, watertightness, consistent orientation (`oriented()`), volume by the
  divergence theorem; `TrimmedSurface::integrate`, `flipped`, `is_untrimmed`; STEP export via
  gmsh/OCC (`export_step`, `add_shell`); LLM tool `sew`.
- Planar sections (`analysis/section.hpp`): `Plane`, `planar_section`, `planar_sections`,
  `intersect(curve, plane)`; LLM tool `section`.
- Trimmed surfaces (`topology/trimmed_surface.hpp`): loops in the parameter domain, `contains`,
  area by Green's theorem, 3D boundary curves, `split` along surface intersections (chords and
  nested loops); exact gmsh export of trimmed surfaces (`add_trimmed_surface`); LLM tools
  `trim_circle`, `split`.
- Public globally adaptive integrator `integrate` (bounded interval budget, absolute target).
- Surface–surface intersection (`analysis/surface_intersect.hpp`): subdivision seeds with normal-cone loop test,
  predictor–corrector marching with 3D and parameter-space step control, exact knot-line and
  boundary points, verified fit of 3D curves and p-curves, seam/pole-aware de-duplication.
- Banded LU solver for interpolation and least squares (`fitting`), dense fallback.
- Curve–curve and curve–surface intersection: overlapping (coincident) portions detected by local
  projections instead of deep subdivision.
- Curve–curve and curve–surface intersection (`analysis/intersect.hpp`): subdivision with
  tangent/normal-cone leaf tests, box-constrained Gauss–Newton, merging of tangential contact
  clusters; LLM tool `intersect`.
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
- Python: any sequence of three numbers converts to `Vec3`, including NumPy rows
  (`np.array([x, y, z])`, rows of an (n, 3) array) – previously only tuples and lists.
- `coons()` also accepts the two side curves in either order (the docstring promised automatic
  orientation, but swapped side curves raised "corner gap").
- Sphinx: `|...|` in two docstrings (`curvature`, `offset_curve`) was read as an RST
  substitution.
- Trimmed boundary curves were fitted across C1 knot lines of the surface (thousands of points,
  slow sewing); they are now split at the knot-line crossings.
- Adaptive quadrature (arc length, area, `integrate`) could explode on noise-dominated
  integrands (local recursive bisection to depth 50); now globally adaptive with a budget.
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
