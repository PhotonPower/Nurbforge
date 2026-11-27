# Changelog

All notable changes are documented here. Format: [Keep a Changelog](https://keepachangelog.com),
versioning: [SemVer](https://semver.org).

## [Unreleased]

### Added
- Flattening of doubly curved surfaces (#17): `flatten(surface, options)` lays a grid of surface
  points out in the plane so that the geodesic distances between neighbours are kept as well
  as possible (stress majorisation, then Gauss-Newton), and returns the flat pattern on the
  parameter domain of the surface with the strains that remain: principal stretches per grid
  cell, their extremes, the rms strain and both areas. Planes and developable surfaces come
  out without strain; `unroll` remains the exact way for those. LLM tool `flatten`.

## [0.17.0] – 2026-11-27
The remaining surface tools of #17 (two-rail sweeps with several profiles, variable offset,
chamfers between surfaces, unrolling), silhouettes, hidden-line removal, Make2D and outlines
(#95), and analysis of curves and surfaces. The API is not yet stable (0.x).

### Added
- Two-rail sweep with several profiles (#17): `sweep2(profiles, rail1, rail2, params, scaling)`
  blends profiles given where they span the rails, as `sweep` does along one trajectory.
  `Sweep2Scaling.NonUniform` stretches the sections along the rail-to-rail direction only
  (heights are kept); `Uniform` is the previous behaviour.
- Variable surface offset (#17): `offset_surface_variable(surface, distance)` moves a surface
  along its normal by a distance varying over it, given as a function `distance(u, v)` or at
  the four corners of the parameter domain (bilinear). Same fit and checks as `offset_surface`.
  LLM tool `offset` takes `corner_distances`.
- Chamfer between surfaces (#17): `chamfer_surfaces(s1, s2, distance1, distance2, options)`
  returns the ruled surface between the curves at the given distances from the intersection
  curve, measured within the surfaces (geodesics normal to the curve), with the contact curves
  and their p-curves; sides as for `fillet_surfaces`. LLM tool `chamfer` takes two surfaces.
- Developable surfaces (#17): `developability(surface)` reports whether a surface is ruled
  along a parameter direction and how far the normal turns along its rulings;
  `unroll(surface, tolerance)` returns the flat pattern of a developable surface (cylinders,
  cones, general ruled developables) in the plane, without stretching and on the same parameter
  domain, so curves on the surface map to the pattern through their parameters.
- Curve analysis (`analysis/curve_analysis.hpp`), the counterparts of the analysis commands of
  an interactive modeller: `curve_curvature(c, u)` gives the curvature with radius, principal
  normal and centre of the osculating circle; `curvature_graph(c, scale, samples_per_span)`
  the hairs and graph line of a curvature comb, with both values at knots where the curvature
  jumps; `curve_direction(c, plane_normal)` the end points and tangents and, for closed curves,
  the sense of rotation and the enclosed area; `is_closed`, `is_planar`;
  `continuity(a, b)` / `continuity(a, end_a, b, end_b)` the geometric (G0/G1/G2) and parametric
  continuity of a joint with gap, tangent angle and curvatures, and `continuity_breaks(c)` the
  knots where a curve is not curvature continuous; `find_duplicate_curves(curves, tolerance)`
  the groups of geometrically identical curves, whatever their direction or parametrisation.
  Length (`arc_length`) and degree were already available.
- Surface analysis (`analysis/surface_analysis.hpp`): `curvature_map(surface | face)` samples
  the Gaussian, mean and principal curvatures on a grid (knot lines included) and reports
  their ranges, where they are attained and the smallest radius of curvature;
  `draft_angle` / `draft_angle_map(surface | face, pull_direction)` the draft angles for a pull
  direction (negative: undercut); `surface_direction` the parameter directions, the normal,
  closed directions, poles and planarity; `continuity(a, side, b)` the G0/G1/G2 continuity of
  two surfaces along a side, and `edge_continuity(shell)` the same across every edge of a shell
  (smooth edges such as seams and fillet borders); `find_duplicate_surfaces(surfaces,
  tolerance)` the groups of identical surfaces. In Python the maps are NumPy arrays.
- Silhouette curves (#95): `silhouette(surface | face | shell, view)` gives the curves where
  the normal is perpendicular to the view ray, for parallel views (`View::parallel(direction)`)
  and perspective views (`View::perspective(eye)`), as 3D curves with p-curves fitted to a
  tolerance. Clipped to the trimmed region of faces; closed where closed; straight to round-off
  for cylinders and cones seen from the side. Reported with them: the parameters where the
  projection has a cusp, and creases of a surface between a part facing the viewer and one
  facing away. Python: `silhouette(obj, direction=..., eye=...)`.
- Hidden-line removal (#95): `visibility(curves, view, occluders)` splits curves where their
  visibility against shells changes and classifies the pieces. Curves lying on a face (its
  edges, silhouettes) are not hidden by that face; rays through an edge of a shell count as
  hits (the hidden edges of a cube seen along its diagonal).
- Make2D (#95): `make2d(shells, direction)` returns the visible and the hidden lines of shells
  for a parallel view, projected exactly onto a plane normal to the view: edges and
  silhouettes, split by visibility, with duplicates removed (visible lines take precedence)
  and smooth edges such as seams left out (`tangent_edges`).
  With `View::perspective(eye)` (Python `eye=`) the drawing is the central projection onto
  the plane through `origin` normal to the line to the eye, an exact projective map of the
  curves; `outline` takes an eye point as well.
- Visibility events (#95): the curves are cut exactly where their projections cross the
  projected edges and silhouettes of the occluders, so pieces end at these crossings to the
  accuracy of the curve intersection and short hidden stretches are found whatever the
  sampling. The sampling remains for curves passing through faces
  (`VisibilityOptions::events`, `samples`).
- LLM tools `silhouette`, `make2d`, `outline`, `unroll`, `analyze_curve`, `analyze_surface`,
  `continuity` and `find_duplicates`.
- Outline (#95): `outline(shells, direction)` gives the contour of shells seen along a
  direction as closed loops (a `Region`: outer contours counter-clockwise, holes clockwise),
  the boundary of the union of their projections.

### Fixed
- Build: configuring an existing build directory again failed with "Unknown CMake command
  pybind11_add_module" when pybind11 had been fetched (no installed pybind11).

## [0.16.0] – 2026-11-19
Blend and match surfaces with G0/G1/G2 continuity, and fitting of surfaces to scattered points
(#17). The API is not yet stable (0.x).

### Added
- Blend and match surfaces (#17): `blend_surfaces(s0, side0, s1, side1, continuity, tension)`
  builds a G0 (ruled), G1 (cubic) or G2 (quintic) transition surface between two surface sides;
  `match_surface(surface, side, target, target_side, continuity)` modifies a surface near one
  side so that it joins the target with G0/G1/G2 and keeps the opposite side and its cross
  derivatives. Both continue the surfaces as exact reparametrisations across the edges; exact
  NURBS when the weights do not vary across the sides (non-rational surfaces, extrusions,
  revolutions at profile ends), otherwise the boundary data are interpolated to a tolerance.
  LLM tools `blend_surfaces` and `match_surface`.
- Fitting to scattered points (#17): `fit_surface(points, ...)` fits a B-spline surface to
  unordered points by least squares with a thin-plate smoothing term (relative weight), solved
  as a banded system, with parameter correction. Parameters from the best-fit plane (height
  fields), a base surface, or four boundary curves, which are then interpolated exactly.
  Returns `FittedSurface` (surface, max / rms error, parameters). LLM tool `fit_surface`.

## [0.15.0] – 2026-11-19
Vertex blends and variable-radius fillets and chamfers of solid edges (#18), and surface tools
(#17): sweeps with several profiles, rail revolve, tapered extrusion, ribbons, cones to a point,
three-sided patches, exact Coons patches from arcs, rebuild of curves and surfaces. The API is
not yet stable (0.x).

### Added
- Vertex blends (#18): `fillet_edges` / `chamfer_edges` handle chains of different face pairs
  meeting at a vertex of three faces. Two chains with a smooth third edge (tangent chains such
  as the rim of a rounded rectangle) continue each other at their common cross-section; two
  chains at a sharp third edge are mitred along their exact intersection curve (a blend running
  past the sharp edge is trimmed by the face there); three chains get a corner patch: for
  fillets the exact spherical triangle of the ball touching all three faces, for chamfers the
  planar triangle of the contact crossings. All groups are now built at once on the input
  shell and sewn once.
- Variable-size edge blends: `sizes=[(point, size), ...]` (C++ `EdgeBlendOptions::sizes`)
  gives the radius or chamfer distance at points on the edges, blended with a C2 quintic step
  in arc length; cross-sections stay exact arcs of the local radius. `fillet_surfaces` takes a
  `radius_at` function (C++ `SurfaceFilletOptions::radius_at`).
- LLM tool `round_edges`: `'*'` wildcards for the face labels (e.g. every edge of the top face)
  and optional variable `sizes`.
- Surface tools (#17): `sweep(profiles, trajectory, params)` (several profiles, blended
  between their positions), `rail_revolve(profile, rail, axis_point, axis_dir, scale_height)`,
  `extrude_tapered(profile, direction, angle)` (draft angle, mitred corners), `ribbon(curve,
  width, plane_normal)`, `extrude_to_point(curve, apex)` (exact cones), `triangle_patch(c0, c1,
  c2)`, `rebuild_curve` / `rebuild_surface` (prescribed control net, reported max error). LLM
  tools `rail_revolve`, `extrude_tapered`, `ribbon`, `extrude_to_point`, `triangle_patch`,
  `rebuild`; `sweep` takes `profiles`.
- `coons` is exact for rational boundaries whose corner weights can be made to agree (e.g.
  arcs): the Boolean sum is taken in homogeneous coordinates.

### Changed
- Edge blends extend polynomial face surfaces linearly (tangent planes) beyond their
  boundaries instead of naturally; natural extensions of fitted surfaces with short spans bent
  away and broke blends near chain ends. Rational surfaces (conics) are still extended
  naturally.
- Edge blend errors at a vertex name its position; a chain that continues tangentially into an
  edge that is not selected is reported as such.

## [0.14.0] – 2026-11-19
Curve tools (chamfers, best-fit lines and planes, projection onto planes, tween curves, tangent
arcs and lines, curves on surfaces with geodesic offsets) and fillets and chamfers of solid
edges (single edges, chains and closed loops, convex and concave). The API is not yet stable
(0.x).

### Added
- Curve tools (#16): `chamfer_curves(a, b, distance_a, distance_b, trim, corner)` and
  `chamfer_corners(curve, distance)` (straight cuts, counterpart of the fillets; the corner
  handling of `fillet_corners` is shared); `fit_line` / `fit_plane` (orthogonal regression,
  with max and rms deviation); `project_to_plane(curve | surface, plane, direction)` (exact
  parallel projection); `tween_curves(a, b, count)`; `make_arc_tangent(start, tangent, end)`;
  `tangent_lines(point, curve)` and `tangent_lines(a, b)` (all common tangents of two planar
  curves); curves on surfaces: `interpolate_on_surface(surface, points, degree, uv)` and
  `offset_on_surface(surface, curve, distance)` (geodesic offset, RK4 on the geodesic
  equations). LLM tools `chamfer`, `chamfer_corners`, `tween`, `project_to_plane`,
  `offset_on_surface`.
- Edge blends on solids (#18): `fillet_edges(shell, edges, radius)` and
  `chamfer_edges(shell, edges, distance)` round or chamfer edges, edge chains and closed edge
  loops (cylinder rims, bores) of a closed shell; convex edges lose material, concave edges gain
  it. The adjacent faces are cut back along the exact contact curves, the end faces of open
  chains are rerouted along the exact end curves; new faces are labelled `fillet` / `chamfer`.
  `edges_between(shell, label_a, label_b)` selects edges by face labels. LLM tool
  `round_edges`. `fillet_surfaces` takes `near` to choose among several spines.
- C++: `Plane` moved to its own header `nurbforge/geometry/plane.hpp` (still included by
  `analysis/section.hpp`).

## [0.12.1] – 2026-11-18
Bug fix: critical heights perpendicular to a cylinder axis (#83).

### Fixed
- `critical_heights` lost height extrema that fall exactly on a sample point of its sign-change
  search, where the derivative is exactly zero (#83): since 0.7.0 the lid edges of a cylinder are
  exact rational circles whose knots lie on the extrema along x and y, so the heights of a cylinder
  perpendicular to its axis came back as one value (FourierGrating's CAD slicer: bounding boxes
  collapsed to a line). An exact zero now counts when the nearest nonzero samples on both sides
  have opposite signs; a run of zeros (a straight edge normal to the direction) counts once. Same
  for edges and for faces whose height depends on one parameter.

## [0.12.0] – 2026-11-18
Solids for simulation (split by a plane with labelled cut faces, cap, thicken, hollow, pipes,
defeaturing), parallel parameter studies from Python (expensive operations release the GIL,
all geometry types picklable), 2D profiles (planar region Booleans, blend curves, projection
onto surfaces), a user manual (PDF) with a generated command reference, a compact API
reference for language models and coding agents shipped in the wheel, and fixes for
Booleans through poles and seams, `split` with crossing cuts, `closest_point` at kinks,
curvature at sphere poles and the installed CMake package. The API is not yet stable (0.x).

### Added
- Compact API reference for language models and coding agents (`REFERENCE.md`, ~20k tokens:
  conventions, a worked example, every public function, class and LLM tool with signature and
  description), shipped in the wheel: `python -m nurbforge` prints it, `nurbforge.reference()`
  returns it; attached to the GitHub release as `nurbforge-reference.md`, indexed by `llms.txt`.
  Generated by `docs/manual/gen_agent_reference.py` (`--check` in CI and the tests).
- 2D profiles: `region_boolean(a, b, op)` (union / intersection / difference of planar regions
  with exact curves; regions are loops of curves, nested loops are holes, `PlanarSection`s are
  accepted) and `region_area`; `blend_curves(a, end0, b, end1, continuity, tension)` (G1 cubic or
  G2 quintic transition); `project_curve(curve, surface, direction, tolerance)` returning the 3D
  curve and its p-curve. LLM tools `blend` and `project`.
- Parallel parameter studies in Python: the expensive operations (Booleans, plane splits, all
  solid operations, intersections, closest points, sections, offsets, sweeps, lofts, Coons and
  Gordon patches, optical surface fits, lengths, areas, volumes, mass properties, sewing) release
  the GIL, so they run in parallel threads; all geometry types (points, transforms, curves,
  surfaces, trim loops, trimmed faces, shells, planes, tolerances, options, optical sags) can be
  pickled for process pools and `copy.deepcopy`. `Shell.tolerance()` (C++) gives the sewing
  tolerance. A ThreadSanitizer test runs constructive operations concurrently.
- `split(shell, plane, tolerance=1e-6, cut_label="cut")` (C++ and Python): splits a closed solid
  by a plane into the closed parts below and above; the new faces in the plane are labelled
  (e.g. for symmetry boundary conditions).
- `cap(shell, tolerance=1e-6, label="cap")`: closes the planar openings of a shell with planar
  faces trimmed exactly by the opening's edges (nested openings become holes).
- `thicken(shell, thickness, symmetric=False, tolerance=1e-6, wall_label="wall")`: gives an
  open, smooth shell a wall thickness (faces offset with their trim loops, ruled side walls).
- `make_pipe(path, radius, end_radius=-1)` / `make_pipe(path, radius_function)`: solid pipe
  along a path (rotation-minimising frame, constant, tapered or function radius), faces labelled
  `wall`, `inlet`, `outlet`; `pipe_surface(path, radius_function, options)` for the wall alone.
- `remove_holes(shell, max_size)`: defeaturing - fills through and blind holes smaller than
  `max_size` and deletes the faces of the removed features.
- `hollow(shell, thickness, open_labels=[])`: hollows a smooth solid to a wall (inner void, or
  open at the faces with the given labels).
- `merge_faces(shell)`: defeaturing - merges adjacent faces with the same label on the same
  surface (pieces of a split face, coplanar planar faces) into one exactly trimmed face.
- LLM tools for the solid operations: `split_solid`, `cap`, `thicken`, `hollow`, `pipe`,
  `remove_holes`, `merge_faces`.
- User manual (LaTeX/PDF, `docs/manual`): installation, concepts and every part of the API in
  18 chapters with 40 tested example scripts, plus a command reference generated from the
  installed package (every class, method, property, enumeration and function of `nurbforge`,
  `nurbforge.mesh` and `nurbforge.llm`, and all LLM tools). `python docs/manual/build.py` builds
  it; CI builds it for documentation changes and attaches `nurbforge-manual.pdf` to every
  release. A test runs all examples and fails if a public API entry has no description.
- `TrimmedSurface.transformed(t)` (loops and label kept, reflections keep the normal) and
  `transform * obj` for curves, surfaces and trimmed faces in Python (was shells only).
- Python: `split(face, cuts)` for trimmed faces (general planar arrangement: cuts clipped to the
  face, may end on trim loops or other cuts and may cross).

### Fixed
- `surface_curvature` (and therefore `offset_surface`, `thicken`) threw "cannot normalize a
  zero-length vector" near the poles of (flipped) spheres: the principal-direction vector of a
  nearly umbilic point is tiny but valid there; directions are now scaled by their own length.
- Installed CMake package: `NurbforgeConfig.cmake` now finds its `Threads` dependency
  (`find_dependency`), so `find_package(Nurbforge)` works without a preceding
  `find_package(Threads)`, and a version file allows `find_package(Nurbforge 0.9)` (#67). CI
  builds a consumer project (`tests/install`) against the installed package on all platforms.
- Tutorial: NumPy arrays of three numbers (e.g. rows of an `(n, 3)` array) are accepted as
  points; the tutorial said they were not (#70).
- `closest_point` on curves and surfaces with kinks (C0 knots, e.g. polylines, rectangle
  outlines and their extrusions) could return a point on the wrong side of a corner: Newton ran
  on the whole curve and took the one-sided derivative of the next span at the segment end. It
  now runs on the Bézier segment / patch itself. This also fixes `hatch(contours=...)` for
  loops made of single kinked curves, which failed with "contour ... does not close" (#69).
- `split(surface, pcurves)` with curves that cross or touch each other (e.g. a cross, or
  chords meeting at a corner) returned the unsplit surface; such cuts now go through the general
  planar arrangement (#68).
- Booleans with cut planes through the poles or along the seam of a closed surface (a sphere
  minus a box with its corner at the centre; a cylinder cut along its seam line) were open with
  GCC and MSVC and closed with Apple Clang (#63). An intersection curve along a seam has its
  p-curve on the domain boundary, where the in/out test depended on the sign of rounding errors
  (FMA contraction); points on a seam now count as inside if the face covers either side of it,
  and intersection curves that repeat another one in 3D (one per side of the seam) are used once.
  A test covers all three operations; the whole C++ suite also passes with `-ffp-contract=off`.

### Changed
- All C++ sources are formatted with clang-format 23.1.5; CI now enforces the format (it was
  report-only).

## [0.9.0] – 2026-11-17
Meshing for solvers: local mesh sizes by face label, curvature-based sizing and BEM-ready
surface meshes (watertight, outward-oriented, as numpy arrays), plus a fix for the gmsh
transfer of Boolean results with trimmed cylinders, cones, spheres and tori. Wheels for
CPython 3.14. The API is not yet stable (0.x).

### Added
- BEM-ready surface meshes: `nurbforge.mesh.bem_mesh(shell, size, order, curvature, refine, path)`
  and `GmshSession.surface_mesh()` return a `SurfaceMesh` (numpy nodes, linear or quadratic
  triangles, groups per face label, normals, areas, volume) that is checked to be watertight
  and oriented outward – into the void for cavities. The orientation is also applied to the
  gmsh model, so written `.msh`/`.stl` files carry it; violations raise `MeshTopologyError`.
- Local mesh sizes: `GmshSession.refine(labels_or_tags, size, distance=0)` – exact size on the
  faces (by face label or gmsh surface tag) and, with `distance`, linear growth to the global size
  into the volume; `mesh(curvature=n)` for n elements per 2π of curvature.
- `GmshSession.add_shell` (remembers face labels) and `GmshSession.face_tags()`: current gmsh
  surface tags per face label, robust against renumbering by duplicate removal.
- Wheels for CPython 3.14 (all platforms).

### Changed
- Releases are published on PyPI (`pip install nurbforge`) via trusted publishing from the
  wheels workflow, in addition to the GitHub release.
- Wheels workflow: `actions/upload-artifact` v4 → v7 (runs on Node.js 24; Node.js 20 is
  deprecated on GitHub runners).

### Fixed
- `nurbforge.mesh.add_shell` transferred Boolean results with trimmed closed faces (bores and
  slots in cylinders, cut spheres, cones, tori) incorrectly to gmsh: faces vanished or the volume
  was wrong, so meshes of such parts were invalid. Shells are now handed over through the native
  STEP writer and OCC's STEP reader, which keeps shared and seam edges; `add_shell` returns the
  surface tags of this shell only, in face order.

## [0.8.0] – 2026-11-14
macOS as a third supported platform (Apple Silicon and Intel, macOS ≥ 14, with binary wheels),
face labels that survive Booleans, transforms and STEP, and gmsh physical groups from them. The
API is not yet stable (0.x).

### Added
- Face labels (#18, stable identification across operations): `TrimmedSurface::label()` /
  `with_label()`, `Shell::with_face_labels()` / `faces_labeled()`; kept by `flipped`,
  `with_hole`, `split`, Booleans (each piece of a split face inherits its label), `oriented`,
  `transformed`, and written to / read from STEP as the `ADVANCED_FACE` name. The solid
  primitives label their faces (box `bottom`/`top`/`front`/`back`/`left`/`right`, cylinder and
  cone `side`/`bottom`/`top`, sphere and torus `surface`).
- `nurbforge.mesh.add_shell(shell, physical_groups=True)`: one named gmsh physical surface group
  per face label and a physical volume, for boundary conditions after Booleans.
- Wheels workflow: a `release` job publishes the GitHub release for every `v*` tag (wheels,
  sdist, the tag's CHANGELOG section as notes); existing tags can be released by running the
  workflow manually with `release_tag`.
- macOS support (Apple Silicon and Intel, macOS ≥ 14): CI job with Apple Clang (`-Werror`, C++
  and Python tests), wheels for `arm64` and `x86_64`, build guide `docs/build/macos.md`.

### Fixed
- STEP reader on macOS: numbers are parsed with `strtod_l` in the "C" locale instead of
  `std::from_chars(double)`, which Apple's libc++ only provides from macOS 26. A test checks
  that a decimal-comma `LC_NUMERIC` (e.g. `de_DE`) does not affect parsing.
- The gmsh physical-group test is skipped instead of failing when gmsh cannot load its system
  libraries (manylinux wheel tests).

## [0.7.0] – 2026-11-10
Writing paths for two-photon lithography (hatching of section loops), solid primitives with mass
properties, curve constructors, and much faster Booleans of STEP shells. The API is not yet
stable (0.x).

### Added
- `hatch(loops, direction, spacing, options)`: hatch and contour lines of the region bounded by
  section loops, e.g. the writing paths of a layer in two-photon lithography (#44). Hatch lines
  n · p = phase + k · spacing are intersected exactly with the loop curves (sign changes of the
  rational height of every Bézier piece, transversal roots refined to machine precision); the
  region follows the loop orientation (winding number), so holes and islands in holes are right.
  Tangential contacts and edges lying on a line do not split it (no zero-length segments).
  Segments carry their line index, the outer loop of their region and the loops they start and
  end on. Contours are the boundaries of the region shrunk by each distance: the raw offsets of
  all loops (arcs at reflex corners, parts beyond the radius of curvature left out) are split at
  all intersections, and the fragments at the full distance from the loops and inside the region
  are chained – so contours are trimmed at convex corners, split at necks narrower than twice the
  distance, merge where holes come close (`HatchContour.loops`) and vanish where the region is too
  thin. Checked against shapely's `buffer(-d)` on ~2000 random polygons and spline loops. Python:
  `nf.hatch(loops, direction, spacing, phase=…, contours=[…], min_length=…)` →
  `(segments, contours)`; example `examples/python/11_slice_and_hatch.py`.
- Solid primitives as closed shells with exact faces and outward normals: `make_box`,
  `make_sphere`, `make_cylinder`, `make_cone` (frustum, or cone with apex), `make_torus` (#18).
- `Shell::mass_properties()`: volume, centroid, inertia tensor about the centroid (unit density),
  principal moments and axes – face integrals by the divergence theorem relative to the bounding
  box centre, voids subtracted (#18).
- `Shell::transformed(transform)` (Python also `transform * shell`): affine maps of shells; mirrors
  flip the faces so normals stay outward. The LLM tools `translate`, `rotate`, `mirror` now work
  on shells, and the new tools `solid` and `mass_properties` create and measure solids.
- `TrimmedSurface.flipped()` in Python.
- Curve constructors (#16): `make_circle_3pt`, `make_arc_3pt` (exact), `make_polygon` (regular
  n-gon), `make_rectangle_curve` (closed outline), `make_helix` (cylindrical or conical helix,
  planar Archimedean spiral, left- or right-handed; parameter = number of turns; C2 quintic
  Hermite spline verified against the exact curve, span count from the h^6 convergence).
- `fillet_corners(curve, radius)`: rounds every corner of a planar curve with exact arcs
  (`fillet_curves` per corner, also the closing corner) and joins it into one G1 curve.
- LLM tools `make_polygon`, `make_helix`, `round_corners`.

### Fixed
- Boolean difference returned open shells for parts with coplanar lids (#45):
  - `split` (face arrangement) probed each region only against its own boundary cycle; for a cut
    around an existing hole (re-drilling a ring with a larger hole) the probe of the outer region
    landed in the old hole, the region was taken for outside the face and the remaining annulus
    was lost. The probe ray now stops at the first crossing of any cycle.
  - On a lid shared with the other solid the cut circle arrived twice – as a fragmented
    intersection curve (it lies on the boundary of the other solid's side face) and as the
    projected loop of the coincident lid. `prune_cuts` kept the fragments and dropped the whole
    loop, so the cut stayed open and the lid was not split (and the arrangement spent ~20 s on the
    overlapping curves). Longer cuts are now kept first. Coaxial cylinders with both lids
    coplanar: open shell in 25 s → closed, exact volume, 2.8 s.
- `TrimmedSurface::boundary_curves` merged break points only closer than 1e-12, but `join` needs
  the joint at least the parametric resolution (1e-10) inside the joined domain: a trim curve
  with a break ~1e-11 from a knot threw "remove_knot: … is not an interior knot" (hit by the
  coaxial overlapping cylinders after the fix above). Breaks closer than 10 · 1e-10 are merged,
  the domain ends kept exactly.
- Booleans of shells read from STEP were ~10x slower than of the originals (#46): every edge
  went through cubic interpolation twice and came back with hundreds of control points.
  - `TrimmedSurface::boundary_curves` (used by `sew` and the STEP writer) is exact where it can
    be: on a planar parallelogram the edge is the affine image of the trim curve (a trim circle
    stays a rational circle), and a straight trim segment along an iso line gives the iso-curve.
    Both keep the trim curve's parametrisation. A cylinder with lids: STEP file 155 kB → 6 kB.
  - STEP reader: edges lying on an iso line of their face (circles at the ends of a cylinder,
    seams, meridians) become straight trim segments instead of interpolants (the end circle on a
    cylinder's side face: 1521 control points → 2).
  - `projected_loops` (cuts of coincident faces) interpolates piecewise between the C1 knots of
    the trim curve; a single cubic across them ran into the sample limit (2049 control points
    for a circle).

  Issue cases (intersection): disjoint cylinder/ring 2.6 s / 23 s (from STEP) → 0.3 s / 0.2 s;
  nested cylinders with coplanar lids 34 s / 73 s → 0.2 s / 0.2 s. `nf_stress 3`: 1.03 → 0.85 s
  per case.

## [0.6.0] – 2026-11-05
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
- Slicing of B-rep shells (`analysis/shell_section.hpp`, Python `section`, `critical_heights`,
  `is_extrusion_band`): `section(shell, plane)` returns closed loops of section pieces (one per face
  crossing, clipped to the trimmed region, chained within the join tolerance), oriented with the
  solid on the left (outer loops counter-clockwise, voids clockwise), with the signed area, the
  source face, its p-curve and the outward normal at every point. Faces whose height along the
  normal depends on one parameter (extrusions, revolutions about the normal, planes) give exact
  isoparametric curves; others are traced to 1e-9 of the model size. `critical_heights` finds
  vertices, planar faces normal to the direction and the height extrema on faces (Newton) and
  edges; `is_extrusion_band` checks a band for walls parallel to the direction. Tests: box, cylinder,
  sphere, frustum, torus (both axis directions), box with a spherical cavity, tilted traced section,
  concurrent queries (ThreadSanitizer CI job).
- STEP: length unit micrometre (`StepLengthUnit::Micrometre`, `SI_UNIT(.MICRO.,.METRE.)`, Python
  `unit="um"`) for writing and reading; shell names (`StepOptions::shell_names`, Python `names=`)
  written as solid names and read back by `read_step_named` / `from_step_named` (Python: lists of
  `(name, Shell)`; fallback: the product name); `StepReadOptions::file_unit_metres` reports the
  length unit declared in the file.
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
