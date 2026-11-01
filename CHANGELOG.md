# Changelog

All notable changes are documented here. Format: [Keep a Changelog](https://keepachangelog.com),
versioning: [SemVer](https://semver.org).

## [Unreleased]
### Fixed (found by the new Boolean stress test, `benchmarks/stress.cpp`)
- `split` (face arrangement) dropped regions thinner than 1e-4 of the parameter domain (its
  region probe was nudged by a fixed 1e-4): a box face overlapping a cylinder by 1e-5 gave open
  Boolean results. The probe is now the midpoint of an inward ray to the first crossing of the
  cycle. Touching runs of intersection curves along a face boundary are recognised within the
  tolerance instead of 10 tolerances (a 1e-5 circle segment was taken for a touching run).
- Curve/surface intersection of a curve running just off a surface (a line next to a tilted
  plane, the patch boundary next to a line of contact between a cylinder and a rotated box
  face) subdivided down to the size limit along the whole curve and took minutes; the signed
  distance is now scanned and its sign changes refined by bisection.
- Intersection curves touching a domain edge: the two ends at the touching point may overlap
  slightly (their position along the tangent is only determined to ~sqrt(2 r tol)) and are now
  joined as well (cylinder inside a rotated box touching two sides).
- `TrimmedSurface::contains` decided near curved trim loops by the loop's polygon (a trim circle
  deviates from its 64-gon by 6e-4): intersection curves ending on a trim circle were clipped
  short, their cuts did not reach the loop and the face was not split (e.g. a box cutting a
  cylinder lid). Within the chord-deviation band the side is now decided on the exact curve.
- `closest_point` on surfaces split degenerate patches at poles below the parametric resolution
  and threw ("split: parameter … is not inside"); a direction that is too narrow is no longer
  split (also for curves).
- Surface/surface intersection: the pole-jump detection (0.5.x) also fired for ordinary long
  steps on small faces and cut straight intersection lines short; a jump must now exceed the
  predicted parameter change.

### Added
- Booleans with tangential contacts: spheres touching at a point (outside, inside, on the seam,
  at a pole), spheres and cylinders touching box faces along a point or a line, a cylinder
  inside a box touching two of its sides and two cylinders side by side give closed results
  with exact volumes (tests `[tangent]`). The surface/surface intersection reports touching
  contacts as points (micrometre curves around a contact on a seam, fragments along a line of
  contact between curved surfaces) and joins curves that touch a domain edge; the Boolean
  classifies touching faces at an interior point clear of the contact and ignores intersection
  curves that only touch a face boundary.
- Booleans: overlaps thinner than the tolerance count as contacts (the face pair is not cut when
  its intersection curves bound a lens thinner than the tolerance), e.g. two spheres
  overlapping by 1e-7 give two touching lumps and an empty intersection.
- Boolean stress test: six families with exact tangential contact (`nf_stress 20 touch`).
- Native STEP reader: `read_step(path)` / `from_step(text)` (C++ `io/step.hpp`, Python
  `nf.read_step(path, unit=...)`, LLM tool `import_step`) return one `Shell` per
  `MANIFOLD_SOLID_BREP` / `BREP_WITH_VOIDS` / surface-model shell. Part 21 parser (complex
  instances, typed parameters, comments); B-spline curves and surfaces (also rational, Bézier,
  (quasi-)uniform); `PLANE`, `CYLINDRICAL_SURFACE`, `CONICAL_SURFACE`, `SPHERICAL_SURFACE`,
  `TOROIDAL_SURFACE` as exact NURBS patches covering the face; `LINE`, `CIRCLE`, `ELLIPSE`,
  `HYPERBOLA`, `PARABOLA`, `POLYLINE` edges (also via `SURFACE_CURVE` / `SEAM_CURVE` /
  `TRIMMED_CURVE`). Trim loops are rebuilt in each face's parameter space (exact for planes,
  otherwise projection + interpolation verified to tolerance / 10, with seams, poles and loops
  around closed surfaces). Length units (SI prefixes, inch) are converted to the requested unit.
  Round trips of the writer's output reproduce volumes to < 1e-6; OCC-written boxes, cylinders,
  cones, spheres, tori, drilled and filleted boxes match the exact volumes.
- `Shell::sew` rejects curve pairs by control-point bounding boxes before projecting (sewing the
  26 faces of a filleted box: 7.7 s → 0.4 s).
- Boolean stress test `nf_stress` (seven surface families, volume-identity oracle; ctest
  `stress_smoke`).
- Surface/surface intersection: the boundary seeds of the leaf patch pairs (90 % of the time)
  are computed on all hardware threads and merged in sequential order (identical results);
  two spheres 174 → 77 ms, two-sphere Boolean 351 → 240 ms. The core library links
  `Threads::Threads`.
- LLM constraint solving: tool `solve` – a recipe of tool calls with unknowns (`"$r"`,
  `"${2*r}"`, whitelist expression interpreter), bounds and targets on measured properties; a
  bounded Levenberg–Marquardt iteration re-runs the recipe in a scratch workspace and executes
  it for real once the targets are met (`nurbforge.llm.constraints`). The system prompt tells
  the model to use it for indirectly given dimensions.
- Bodies with several shells: `Shell::components()` (Python `Shell.components()`); `oriented()`
  orients closed components by nesting depth, so voids keep inward normals (hollow solids have
  the correct volume and work as Boolean inputs; lumps of a Boolean result are separate
  components).

### Fixed
- `oriented()` of a shell with several components only checked the sign of the total volume,
  so the relative orientation of the components (e.g. solid and void) was arbitrary.
- Surface/surface intersection through the poles of a surface of revolution (e.g. a plane through
  a sphere's axis): boundary points at a pole took an arbitrary parameter, which generated
  spurious knot-line points and made the fitted curve run back and forth (wrong Booleans); pole
  end points now continue the neighbouring parameter and lie exactly on the boundary; marching
  steps against the direction of travel are rejected and a collapsing step next to a pole ends
  on it; corrected midpoints must lie between their neighbours.
- Trim curves running along a knot line of the surface (e.g. such meridians) no longer produce a
  crossing at every rounding-level sign change (zero-length pieces in `boundary_curves`).
- Face splitting: curve ends merged into one vertex are moved onto it, so curves meeting within
  the merge distance (e.g. at a pole) form exactly closed loops.

## [0.5.0] – 2026-11-01
Roadmap milestones 0.2–0.5: NURBS toolbox, surface construction and freeform optics,
intersections and trimming, B-rep shells, Boolean operations and STEP export. The API is not yet
stable (0.x).

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
- Benchmarks: `benchmarks/bench.cpp` (`-DNF_BUILD_BENCHMARKS=ON`, built and smoke-tested in the CI
  preset), reference timings and an OpenMP experiment in `docs/benchmarks.md`.
- LLM verification: tools `inspect` (measured properties and detected issues of an object) and
  `check` (measured property vs. expected value, PASS/FAIL); `GeometryAgent` runs verification
  rounds after the final answer (`verify_rounds`, default 1; CLI `--verify-rounds`): objects
  created or changed for the request are measured, issues such as open shells or degenerate
  geometry are flagged, and the model repairs or confirms. Module `nurbforge.llm.verify`.
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
