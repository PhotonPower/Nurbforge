# Python API overview

```python
import nurbforge as nf
```

## Types

| Class | Key members |
|---|---|
| `Vec3(x, y, z)` | `x, y, z`, `+ - *`, iterable; tuples/lists convert implicitly |
| `NurbsCurve(degree, knots, control_points, weights=[])` | `evaluate(u)` / `c(u)`, `derivatives(u, order)`, `tangent(u)`, `sample(n) -> ndarray (n,3)`, `insert_knot(u, times)`, `refine_knots(X)`, `decompose_bezier()`, `elevate_degree(times)`, `reduce_degree(tolerance)`, `try_remove_knot(u, times, tolerance) -> (curve, removed)`, `remove_knot(u, times, tolerance)`, `split(u) -> (left, right)`, `sub_curve(u0, u1)`, `clamped()`, `is_clamped`, `reversed()`, `transformed(T)`, `domain`, `knots`, `weights`, `control_points`, `is_rational` |
| `NurbsCurve.from_points(degree, pts)` | clamped uniform B-spline |
| `NurbsSurface(du, dv, ku, kv, nu, nv, pts, weights=[])` | `evaluate(u,v)` / `s(u,v)`, `derivatives(u, v, order)`, `first_derivatives`, `normal` (limit normal at poles), `iso_curve_u/v`, `sample(nu,nv) -> ndarray (nu,nv,3)`, `refine_knots_u/v(X)`, `elevate_degree_u/v(times)`, `reduce_degree_u/v(tolerance)`, `remove_knot_u/v(u, times, tolerance)`, `try_remove_knot_u/v`, `split_u/v`, `sub_surface(u0, u1, v0, v1)`, `clamped()`, `decompose_bezier()`, `transformed`, `transposed` |
| `Transform` | `translation`, `rotation(origin, axis, angle)`, `scaling`, `mirror(point, normal)`, `T1 * T2`, `apply`, `apply_vector` |

## Functions

| Group | Functions |
|---|---|
| Primitives | `make_line`, `make_polyline`, `make_arc`, `make_circle`, `make_bilinear`, `make_rectangle` |
| Construction | `extrude`, `revolve`, `ruled`, `reparametrize`, `loft(sections, degree_v=3, method)` / `loft(sections, degree_v, v_params)`, `loft_parameters`, `make_compatible(curves)`, `sweep(profile, trajectory, SweepOptions(frame, twist, tolerance))`, `sweep2(profile, rail1, rail2, options)`, `sweep_translational(profile, trajectory)`, `coons(c0, c1, d0, d1, tolerance)`, `gordon(u_curves, v_curves, tolerance)` |
| Analysis | `closest_point(obj, point, tolerance=Tolerance())` → `CurvePoint(u, point, distance)` / `SurfacePoint(u, v, point, distance)`; `closest_point_local(obj, point, u0[, v0])` |
| Properties | `curvature(curve, u)`, `torsion(curve, u)`, `length(curve)`, `arc_length(curve, u0, u1)`, `parameter_at_length(curve, s)`, `surface_curvature(surface, u, v)` → `SurfaceCurvature(gaussian, mean, k1, k2, normal, dir1, dir2)`, `area(surface[, u0, u1, v0, v1])` |
| Fitting | `interpolate_curve(points, degree=3, method=Parametrization.Centripetal, end_derivatives=None)`, `approximate_curve(points, degree, n_ctrl)`, `approximate_curve_tol(points, degree, tolerance)`, `interpolate_surface(grid (nu,nv,3), degree_u, degree_v)`, `approximate_surface(grid, degree_u, degree_v, ctrl_u, ctrl_v)`, `data_parameters(points, method)` |
| Intersections | `intersect(curve, curve)` → `[CurveCurveIntersection(s, t, point, distance)]`, `intersect(curve, surface)` → `[CurveSurfaceIntersection(t, u, v, point, distance)]`, `intersect(surface, surface, SurfaceIntersectionOptions(tolerance, chord_tolerance, fit_tolerance))` → `SurfaceIntersection(curves=[IntersectionCurve(points, closed, curve, pcurve1, pcurve2)], points)` |
| Sections | `Plane(point, normal)`, `planar_section(surface, plane)` → `SurfaceIntersection` (pcurve2 in in-plane coordinates), `planar_sections(surface, plane, spacing, count)`, `intersect(curve, plane)` |
| Tolerances | `Tolerance(length=1e-10, parametric=1e-10, angular=1e-9)` |
| Extrapolation | `extend(curve, CurveEnd.Start/End, length, Extension.Natural/Linear)`, `extend(surface, SurfaceSide.UMin/UMax/VMin/VMax, length, type)` |
| Offsets | `offset_curve(curve, distance, plane_normal=(0, 0, 1), tolerance=1e-7)` (planar, cubic, parameters kept), `offset_surface(surface, distance, tolerance=1e-7)` (along the normal, quintic); raise `GeometryError` where the distance reaches a radius of curvature |
| Fillets | `fillet_curves(a, b, radius, plane_normal=(0, 0, 1), side_a=0, side_b=0, trim=True, corner=None, tolerance=1e-7)` → `CurveFillet(arc, center, param_a, param_b, point_a, point_b, trimmed_a, trimmed_b)` (exact arc; sides +1 = left of the curve direction seen from `plane_normal`, 0 = automatic corner); `fillet_surfaces(s1, s2, radius, side1=1, side2=1, tolerance=1e-6, offset_tolerance=0, max_sections=1025)` → `SurfaceFillet(surface, spine, contact1, contact2, pcurve1, pcurve2)` (rolling ball; side +1 = ball on the side of S_u × S_v) |
| Optics | `Asphere(curvature, conic, even, odd)`, `asphere_surface(a, r_max, tolerance)`, `asphere_profile`, `ZernikeSag(base, norm_radius, indexing, coeffs)`, `zernike_surface`, `zernike(j, indexing, rho, theta)`, `zernike_nm`, `XYPolynomialSag(base, norm_radius, terms)`, `xy_polynomial_surface`, `sag_surface(f, x0, x1, y0, y1, degree, tolerance)` – see [optics.md](../optics.md) |
| Trimming | `TrimmedSurface(surface[, outer, holes])` -> `contains`, `area`, `boundary_curves`, `with_hole`; `TrimLoop(curves)`, `rectangle_loop`, `circle_loop`, `polygon_loop`; `split(surface, cutter)` / `split(surface, pcurves)` |
| Shells | `Shell.sew(faces, tolerance)` → `faces`, `vertices`, `edges` (`ShellEdge(curve, v0, v1, uses)`), `free_edges()`, `non_manifold_edges()`, `is_closed`, `is_consistently_oriented`, `oriented()`, `area()`, `volume()`; mesh: `add_shell`, `export_step(objects, path)` |
| Booleans | `boolean(a, b, BooleanOp.Union/Intersection/Difference)`, `point_in_shell(shell, point)`, `split(trimmed_surface, cuts)` |
| STEP export | `write_step(path, shell_or_shells, unit="mm", name="nurbforge", author="", organization="", timestamp="", uncertainty=1e-6)` (native AP214 writer: closed shells → solids, open shells → surface models; `unit` in `"mm"`, `"cm"`, `"m"`, `"inch"`), `to_step(shell_or_shells, ...)` → file contents as `str` – see [io.md](../io.md) |
| Integration | `integrate(f, a, b, rel_tol, breaks, abs_tol)` |
| Joining | `join(a, b, tolerance)`, `join_u(a, b, tolerance)`, `join_v(a, b, tolerance)` |
| Arrays | `linear_pattern`, `grid_pattern`, `polar_pattern`, `apply_pattern` |
| Knots | `clamped_uniform_knots`, `unique_knots`, `validate_knots` |
| Vector math | `dot`, `cross`, `norm`, `distance` |

## Sub-packages

- `nurbforge.mesh` – `GmshSession`, `add_curve`, `add_surface`, `mesh_surfaces` ([meshing.md](../meshing.md))
- `nurbforge.llm` – `OllamaClient`, `GeometryAgent`, `Workspace`, CLI `nurbforge-llm` ([llm.md](../llm.md))

## Errors

`nurbforge.GeometryError` (subclass of `ValueError`) for invalid geometry.
