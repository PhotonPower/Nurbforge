# nurbforge 0.16.0 — API reference for agents

Headless NURBS CAD kernel (C++20 core, Python bindings): exact curves and surfaces, trimmed
faces, closed shells (B-rep solids), Booleans, fillets, sweeps, fitting, optical surfaces,
planar slicing/hatching, native STEP read/write, gmsh meshing. `import nurbforge as nf`.
This file is the complete public API. Print it with `python -m nurbforge`;
`nf.reference()` returns it as a string.

## Conventions

- Points/vectors: any 3-sequence (tuple, list, 1-D array) or `nf.Vec3`. Results are `Vec3`
  (`v.x`, `x, y, z = v`, `v.to_tuple()`, arithmetic).
- Angles in **radians** (the LLM tools of `nurbforge.llm` use degrees). No fixed length unit;
  STEP I/O declares/converts units (default mm).
- Every constructed curve/surface has parameter domain [0, 1]; the parameter is not arc length.
  Call objects to evaluate: `curve(u)`, `surface(u, v)`.
- Surface control net is row-major: `index = i * num_v + j`, `i` along u.
- `revolve`: u = circumferential, v = profile. `extrude`/`ruled`: u = profile, v = sweep.
  `loft`: u along the sections, v across them.
- `make_circle(center, normal, r)` runs counter-clockwise about `normal`; its start point is
  chosen by the kernel — use `make_arc`/`make_circle_3pt` to control it.
- Trim loops live in the (u, v) space of the surface: outer counter-clockwise, holes clockwise.
- Solids are closed `Shell`s of `TrimmedSurface` faces with outward normals. Faces carry string
  labels (box: bottom/top/front/back/left/right; cylinder: side/bottom/top; ...). Booleans and
  solid operations keep labels; labels become physical groups when meshing.
- Planar sections: loops with the material on the left (outer CCW, holes CW in plane coords).
- All geometry objects are immutable: methods return new objects. They pickle and deepcopy.
  Heavy functions release the GIL; use threads or processes freely.
- `Tolerance(length=1e-10, parametric=1e-10, angular=1e-9)` is the default `Tolerance()`;
  functions taking `tolerance: float` (Booleans, shells, STEP) default to 1e-6 model units.
- Errors: invalid input or failed geometry raises `nf.GeometryError` (a `ValueError`); a failed
  operation leaves its inputs unchanged.
- Optional extras: `pip install "nurbforge[mesh]"` (gmsh, for `nurbforge.mesh`),
  `"nurbforge[llm]"` (Ollama client, for `nurbforge.llm`).

## Typical workflow

```python
import math
import nurbforge as nf

arc = nf.make_arc((0, 0, 0), (1, 0, 0), (0, 0, 1), 10.0, 0.0, math.pi / 2)
dome = nf.revolve(arc, (0, 0, 0), (0, 0, 1), 2 * math.pi)     # NurbsSurface
print(nf.length(arc), nf.area(dome), dome(0.25, 0.5))

box = nf.make_box((-20, -20, 0), (40, 40, 10))                  # closed Shell
bore = nf.make_cylinder((0, 0, -5), (0, 0, 1), 6.0, 20.0)
part = nf.boolean(box, bore, nf.BooleanOp.Difference)
print(part.is_closed, part.face_labels, part.mass_properties().volume)
print(nf.point_in_shell(part, (10, 10, 5)))

sec = nf.section(part, nf.Plane((0, 0, 5), (0, 0, 1)))          # PlanarSection: loops
lines, contours = nf.hatch(sec, (1, 0, 0), 0.5)                  # hatch segments + contours
left, right = nf.split(part, nf.Plane((0, 0, 0), (1, 0, 0)))    # two closed shells

nf.write_step("part.step", part, unit="mm")
shells = nf.read_step("part.step")                               # list[Shell]
# meshing (needs gmsh): nurbforge.mesh.bem_mesh(part, size=2.0, path="part.msh")
```

Free-form: build curves (`make_*`, `interpolate_curve`, `approximate_curve`), surfaces
(`extrude`, `revolve`, `loft`, `sweep`, `sweep2`, `coons`, `gordon`, `interpolate_surface`),
trim them (`split`, `TrimmedSurface(surface, outer_loop)`), sew faces with
`Shell.sew(faces)`, close planar openings with `cap`, give thickness with `thicken`. Check
results with `shell.is_closed`, `shell.free_edges()`, `shell.volume()`; methods are shown
below as `.name(...)`, properties as `.name: type`.

## API

Signature notation: `Vec3` = any 3-sequence; `array_like` = NumPy array or nested lists;
arguments after `*` are keyword-only; `(1) … (2) …` are overloads.

### Functions of `nurbforge`

**Curve primitives**

- `make_line(p0: Vec3, p1: Vec3) -> NurbsCurve` — Straight segment from `p0` (u = 0) to `p1` (u = 1), degree 1.
- `make_polyline(points: Sequence[Vec3]) -> NurbsCurve` — Degree-1 curve through the points in order (one span per segment, parameters at the cumulative chord length, normalised to [0, 1]). At least two points.
- `make_arc(center: Vec3, x_axis: Vec3, y_axis: Vec3, radius: float, start_angle: float, end_angle: float) -> NurbsCurve` — Exact circular arc (rational quadratic) of the given `radius` around `center`, from angle `start_angle` to `end_angle` (radians), measured from `x_axis` towards `y_axis`. The axes need not be unit vectors but must be perpendicular.
- `make_arc_3pt(p0: Vec3, p1: Vec3, p2: Vec3) -> NurbsCurve` — Circular arc from p0 through p1 to p2.
- `make_circle(center: Vec3, normal: Vec3, radius: float) -> NurbsCurve` — Exact full circle (rational quadratic, 9 control points), counter-clockwise about `normal`; the start point lies in the direction `any_perpendicular(normal)`. Domain [0, 1].
- `make_circle_3pt(p0: Vec3, p1: Vec3, p2: Vec3) -> NurbsCurve` — Full circle through three points (starts at p0, runs through p1 and p2).
- `make_arc_tangent(start: Vec3, tangent: Vec3, end: Vec3) -> NurbsCurve` — Exact circular arc from `start`, leaving it in the direction `tangent`, to `end` (sweep between 0 and 360 degrees). Raises GeometryError if `end` lies on the tangent line.
- `make_polygon(center: Vec3, normal: Vec3, radius: float, n: int, x_axis: Vec3 = (0, 0, 0)) -> NurbsCurve` — Closed regular polygon with n vertices on the circle of the given radius (first vertex along x_axis).
- `make_rectangle_curve(origin: Vec3, edge_u: Vec3, edge_v: Vec3) -> NurbsCurve` — Closed rectangle outline (degree 1).
- `make_helix(axis_point: Vec3, axis: Vec3, start: Vec3, pitch: float, turns: float, end_radius: float = -1.0, left_handed: bool = False, tolerance: float = 1e-07) -> NurbsCurve` — Helix / conical helix / planar spiral starting at `start`, parameter = turns; C2 quintic within `tolerance` of the exact curve. end_radius < 0: cylindrical.

**Surface primitives and construction**

- `make_bilinear(p00: Vec3, p10: Vec3, p01: Vec3, p11: Vec3) -> NurbsSurface` — Bilinear patch (degree 1 in u and v) with the corners `p00` (u = v = 0), `p10` (u = 1, v = 0), `p01` and `p11`.
- `make_rectangle(origin: Vec3, edge_u: Vec3, edge_v: Vec3) -> NurbsSurface` — Planar parallelogram patch `origin + u*edge_u + v*edge_v`, u, v in [0, 1]; its normal is `edge_u x edge_v`.
- `extrude(profile: NurbsCurve, direction: Vec3) -> NurbsSurface` — Exact extrusion of `profile` along `direction`: u follows the profile, v in [0, 1] the extrusion (degree 1).
- `revolve(profile: NurbsCurve, axis_point: Vec3, axis_dir: Vec3, angle: float) -> NurbsSurface` — Exact surface of revolution of `profile` about the axis through `axis_point` along `axis_dir`, by `angle` radians (right-hand rule): u is circumferential (rational quadratic), v follows the profile.
- `ruled(c0: NurbsCurve, c1: NurbsCurve) -> NurbsSurface` — Ruled surface between `c0` (v = 0) and `c1` (v = 1); degrees and knot vectors of the two curves are made compatible first. u follows the curves, v is linear.
- `(1) loft(sections: Sequence[NurbsCurve], degree_v: int = 3, method: Parametrization = Parametrization.Centripetal) -> NurbsSurface`  
  `(2) loft(sections: Sequence[NurbsCurve], degree_v: int, v_params: Sequence[float]) -> NurbsSurface` — (1) Skinned surface through section curves (u along the sections, v across them).
- `loft_parameters(sections: Sequence[NurbsCurve], method: Parametrization = Parametrization.Centripetal) -> list[float]` — The v parameters at which `loft` places the sections for the given method (averaged over the section control points), first 0 and last 1.
- `make_compatible(curves: Sequence[NurbsCurve]) -> list[NurbsCurve]` — Common degree, knot vector and domain [0, 1] for all curves.
- `(1) sweep(profile: NurbsCurve, trajectory: NurbsCurve, options: SweepOptions = SweepOptions()) -> NurbsSurface`  
  `(2) sweep(profiles: Sequence[NurbsCurve], trajectory: NurbsCurve, params: Sequence[float] = [], options: SweepOptions = SweepOptions()) -> NurbsSurface` — (1) Sweep a profile (given at the trajectory start) along a trajectory. (2) Sweep with several profiles, each given where it crosses the trajectory (at params[k]; empty: nearest trajectory point); the shape is blended between them (loft) and placed along the trajectory.
- `(1) sweep2(profile: NurbsCurve, rail1: NurbsCurve, rail2: NurbsCurve, options: SweepOptions = SweepOptions()) -> NurbsSurface`  
  `(2) sweep2(profiles: Sequence[NurbsCurve], rail1: NurbsCurve, rail2: NurbsCurve, params: Sequence[float] = [], scaling: Sweep2Scaling = Sweep2Scaling.Uniform, options: SweepOptions = SweepOptions()) -> NurbsSurface` — (1) Two-rail sweep with uniform scaling. (2) Two-rail sweep with several profiles, each given where it spans the rails (at the rail1 parameter params[k]; empty: nearest point of rail1); the shape is blended between them. scaling: Uniform, or NonUniform (stretched along the rail-to-rail direction only, heights kept).
- `sweep_translational(profile: NurbsCurve, trajectory: NurbsCurve) -> NurbsSurface` — Exact translational sweep S(u, v) = C(u) + T(v) - T(v0).
- `pipe_surface(path: NurbsCurve, radius: Callable[[float], float], options: SweepOptions = SweepOptions()) -> NurbsSurface` — Pipe wall: a circle of radius radius(t) (t in [0, 1] along the path) swept with the rotation-minimising frame; u around the circle, v along the path.
- `rail_revolve(profile: NurbsCurve, rail: NurbsCurve, axis_point: Vec3, axis_dir: Vec3, scale_height: bool = False, options: SweepOptions = SweepOptions()) -> NurbsSurface` — Revolve a profile about an axis while its end nearer to the rail follows the rail (radial scaling about the axis; scale_height: also axial scaling about the axis point). E.g. a reflector with a non-round rim.
- `coons(c0: NurbsCurve, c1: NurbsCurve, d0: NurbsCurve, d1: NurbsCurve, tolerance: float = 1e-07) -> NurbsSurface` — Coons patch: S(u,0)=c0, S(u,1)=c1, S(0,v)=d0, S(1,v)=d1 (orientations fixed automatically). Exact also for rational boundaries (arcs) whose corner weights agree.
- `gordon(u_curves: Sequence[NurbsCurve], v_curves: Sequence[NurbsCurve], tolerance: float = 1e-07) -> NurbsSurface` — Gordon surface through a consistent curve network.
- `triangle_patch(c0: NurbsCurve, c1: NurbsCurve, c2: NurbsCurve, tolerance: float = 1e-07) -> NurbsSurface` — Three-sided patch bounded by three curves forming a loop: a Coons patch with c0 as one side and the opposite corner as a degenerate side (v = 1).
- `blend_surfaces(s0: NurbsSurface, side0: SurfaceSide, s1: NurbsSurface, side1: SurfaceSide, continuity: int = 2, tension: float = 1.0, tolerance: float = 1e-07) -> NurbsSurface` — Blend surface between a side of s0 and a side of s1 continuing both with G0 (ruled), G1 (cubic across) or G2 (quintic across, same normal curvature). u runs along the edge of s0, v from s0 (v = 0) to s1 (v = 1); the edge of s1 is mapped linearly onto it in the orientation pairing the nearer points. The cross speed is tension x the mean edge distance. Exact when the weights do not vary across the sides, otherwise the boundary data are interpolated to tolerance.
- `match_surface(surface: NurbsSurface, side: SurfaceSide, target: NurbsSurface, target_side: SurfaceSide, continuity: int = 2, tolerance: float = 1e-07) -> NurbsSurface` — Modify a surface near one side so that it joins the target's side with G0, G1 or G2 continuity: the first continuity + 1 control rows are replaced (after refinement), the opposite side and its cross derivatives are kept. Same parameter directions and domain.
- `extrude_tapered(profile: NurbsCurve, direction: Vec3, angle: float, tolerance: float = 1e-07) -> list[NurbsSurface]` — Tapered (draft) extrusion of a planar profile normal to `direction` (its length is the height): the top is the profile offset by height*tan(angle), outwards for angle > 0. Corners are mitred; returns one ruled surface per smooth piece of the profile.
- `ribbon(curve: NurbsCurve, width: float, plane_normal: Vec3 = (0, 0, 1), tolerance: float = 1e-07) -> list[NurbsSurface]` — Planar ribbon: ruled surfaces between a curve and its offset by `width` (left of the curve seen from the normal's tip), mitred at corners; one surface per smooth piece.
- `extrude_to_point(curve: NurbsCurve, apex: Vec3) -> NurbsSurface` — Cone over a curve: ruled surface from the curve (v = 0) to the apex (v = 1); exact, also for rational curves.
- `rebuild_curve(curve: NurbsCurve, num_ctrl: int, degree: int = 3) -> RebuiltCurve` — Rebuild a curve as a non-rational B-spline with num_ctrl control points and uniform knots on the same domain (least squares, parametrisation kept); returns RebuiltCurve with the max error.
- `rebuild_surface(surface: NurbsSurface, num_u: int, num_v: int, degree_u: int = 3, degree_v: int = 3) -> RebuiltSurface` — Rebuild a surface with num_u x num_v control points and uniform knots on the same domain; returns RebuiltSurface with the max error.
- `reparametrize(curve: NurbsCurve, a: float, b: float) -> NurbsCurve` — The same curve with its domain mapped linearly to [`a`, `b`] (knots scaled; geometry unchanged).

**Solid primitives**

- `make_box(corner: Vec3, size: Vec3) -> Shell` — Axis-aligned box (closed Shell, outward normals) with a corner and edge lengths.
- `make_sphere(center: Vec3, radius: float) -> Shell` — Closed shell of a sphere: one face (revolved half circle, poles at the axis ends, seam edge), outward normals, face label `surface`.
- `make_cylinder(base: Vec3, axis: Vec3, radius: float, height: float) -> Shell` — Cylinder with the base disc centred at `base` along `axis`.
- `make_cone(base: Vec3, axis: Vec3, r0: float, r1: float, height: float) -> Shell` — Frustum: radius r0 at the base, r1 at the top; one radius may be 0 (apex).
- `make_torus(center: Vec3, axis: Vec3, major_radius: float, minor_radius: float) -> Shell` — Closed shell of a torus with the given `major_radius` (centre to tube axis) and `minor_radius` (tube) about `axis`; one face with seams in u and v, label `surface`. Requires `minor_radius < major_radius`.

**Planar regions and curve tools**

- `region_boolean(a: object, b: object, op: BooleanOp, tolerance: float = 1e-09) -> list[list[NurbsCurve]]` — Union / intersection / difference of planar regions in the xy plane. A region is a PlanarSection, or a list of loops: closed curves, lists of curves running head to tail, or SectionLoops (nested loops: holes). Returns the result as a list of loops (lists of exact curves; outer counter-clockwise, holes clockwise).
- `region_area(region: object) -> float` — Area of a planar region (outer loops minus holes; Green's theorem on the exact curves).
- `blend_curves(a: NurbsCurve, end0: CurveEnd, b: NurbsCurve, end1: CurveEnd, continuity: int = 2, tension: float = 1.0) -> NurbsCurve` — Smooth transition from an end of curve a to an end of curve b: G1 (cubic) or G2 (quintic, same curvature vectors). The blend leaves a and enters b; |C'| at the ends is tension x chord.
- `tween_curves(a: NurbsCurve, b: NurbsCurve, count: int) -> list[NurbsCurve]` — `count` intermediate curves between a and b (made compatible, homogeneous control points interpolated linearly at s = k / (count + 1)). Both should run in the same direction.
- `(1) tangent_lines(a: NurbsCurve, b: NurbsCurve, plane_normal: Vec3 = (0, 0, 1)) -> list[TangentLine]`  
  `(2) tangent_lines(a: Vec3, b: NurbsCurve, plane_normal: Vec3 = (0, 0, 1)) -> list[TangentLine]` — (1) Common tangent lines of two coplanar curves (from a to b), or - with a point as `a` - the lines from the point tangent to curve `b`. Returns a list of TangentLine.
- `project_curve(curve: NurbsCurve, surface: NurbsSurface, direction: Vec3 = (0, 0, 0), tolerance: float = 1e-06) -> ProjectedCurve` — Project a curve onto a surface along its normal (direction zero) or along `direction`; returns ProjectedCurve(curve, pcurve) within `tolerance` of the exact projection.
- `(1) project_to_plane(curve: NurbsCurve, plane: Plane, direction: Vec3 = (0, 0, 0)) -> NurbsCurve`  
  `(2) project_to_plane(surface: NurbsSurface, plane: Plane, direction: Vec3 = (0, 0, 0)) -> NurbsSurface` — (1) Exact parallel projection of a curve onto a plane along `direction` (zero: along the plane normal, i.e. orthogonal projection). (2) Exact parallel projection of a surface onto a plane.
- `interpolate_on_surface(surface: NurbsSurface, points: array_like, degree: int = 3, uv: bool = False, tolerance: float = 1e-06) -> ProjectedCurve` — Curve on a surface through points (n, 3): projected onto the surface (or (u, v, 0) parameter points if uv=True), interpolated in the parameter domain; returns ProjectedCurve(curve, pcurve) with the 3D curve within `tolerance` of S(pcurve(t)).
- `offset_on_surface(surface: NurbsSurface, curve: object, distance: float, tolerance: float = 1e-06) -> ProjectedCurve` — Geodesic offset of a curve within a surface: `curve` is its p-curve (u, v, 0) or a ProjectedCurve; positive distances go to the left of the curve seen from the normal S_u x S_v. Returns ProjectedCurve(curve, pcurve) within `tolerance`.

**Joining, extension, offsets, fillets and chamfers**

- `join(a: NurbsCurve, b: NurbsCurve, tolerance: float = 1e-10) -> NurbsCurve` — Join curve b to the end of curve a (gap <= tolerance).
- `join_u(a: NurbsSurface, b: NurbsSurface, tolerance: float = 1e-10) -> NurbsSurface` — Join two surfaces along u: the u = max boundary of `a` must coincide with the u = min boundary of `b` within `tolerance` (same v knots after making them compatible). Returns one surface.
- `join_v(a: NurbsSurface, b: NurbsSurface, tolerance: float = 1e-10) -> NurbsSurface` — Join two surfaces along v (as `join_u` with the roles of u and v exchanged).
- `(1) extend(curve: NurbsCurve, end: CurveEnd, length: float, type: Extension = Extension.Natural) -> NurbsCurve`  
  `(2) extend(surface: NurbsSurface, side: SurfaceSide, length: float, type: Extension = Extension.Natural) -> NurbsSurface` — (1) Extend a curve at one end by an arc length; original parameters are kept. (2) Extend a surface across one side; original parameters are kept.
- `offset_curve(curve: NurbsCurve, distance: float, plane_normal: Vec3 = (0, 0, 1), tolerance: float = 1e-07) -> NurbsCurve` — Planar offset c(t) + d * normalize(plane_normal x c'(t)) as a cubic B-spline fitted to the tolerance (parameters kept). Raises GeometryError at cusps (`|d|` >= radius of curvature on the concave side), corners or non-planar curves.
- `offset_surface(surface: NurbsSurface, distance: float, tolerance: float = 1e-07) -> NurbsSurface` — Offset surface S(u,v) + d * N(u,v) as a quintic tensor-product spline fitted to the tolerance (parameters kept). Raises GeometryError if d reaches a radius of curvature.
- `(1) offset_surface_variable(surface: NurbsSurface, corner_distances: Sequence[float] (4), tolerance: float = 1e-07) -> NurbsSurface`  
  `(2) offset_surface_variable(surface: NurbsSurface, distance: Callable[[float, float], float], tolerance: float = 1e-07) -> NurbsSurface` — (1) Variable offset with the distance given at the corners (u0,v0), (u1,v0), (u0,v1), (u1,v1) of the parameter domain, interpolated bilinearly. (2) Variable offset S(u,v) + d(u,v) * N(u,v) with a smooth distance function of the surface parameters (e.g. a layer with a thickness profile); fit, parametrisation and errors as offset_surface.
- `fillet_curves(a: NurbsCurve, b: NurbsCurve, radius: float, plane_normal: Vec3 = (0, 0, 1), side_a: int = 0, side_b: int = 0, trim: bool = True, corner: Vec3 | None = None, tolerance: float = 1e-07) -> CurveFillet` — Circular fillet of the given radius between two coplanar intersecting curves. side_a / side_b: +1 = centre to the left of the curve's direction seen from plane_normal, -1 = right, 0 = automatic (the corner between the kept parts, i.e. the longer sides of the intersection). Returns a CurveFillet (exact arc, tangent parameters, trimmed curves).
- `fillet_corners(curve: NurbsCurve, radius: float, tolerance: float = 1e-07) -> NurbsCurve` — Round every corner of a planar curve (polyline, polygon, ...) with exact arcs; one G1 curve.
- `fillet_surfaces(s1: NurbsSurface, s2: NurbsSurface, radius: float, side1: int = 1, side2: int = 1, tolerance: float = 1e-06, offset_tolerance: float = 0.0, max_sections: int = 1025, near: Vec3 | None = None) -> SurfaceFillet` — Rolling-ball fillet of radius r between two surfaces. side1 / side2 = +1 put the ball on the side of the normal S_u x S_v of s1 / s2, -1 on the opposite side; `near` selects the spine nearest to a point if there are several. Returns a SurfaceFillet (surface, spine, contact curves, p-curves on s1 / s2).
- `chamfer_curves(a: NurbsCurve, b: NurbsCurve, distance_a: float, distance_b: float | None = None, trim: bool = True, corner: Vec3 | None = None) -> CurveChamfer` — Chamfer between two intersecting curves: a straight cut from the point of a at distance `distance_a` from the corner to the point of b at `distance_b` (default: distance_a), on the kept parts (the longer sides of the intersection, as in fillet_curves). The curves need not be planar. Returns a CurveChamfer (line, parameters, trimmed curves).
- `chamfer_corners(curve: NurbsCurve, distance: float, tolerance: float = 1e-07) -> NurbsCurve` — Chamfer every corner of a curve (polyline, polygon, ...) with straight cuts at `distance` from the corner on both sides; one curve.
- `chamfer_surfaces(s1: NurbsSurface, s2: NurbsSurface, distance1: float, distance2: float | None = None, side1: int = 1, side2: int = 1, tolerance: float = 1e-06, near: Vec3 | None = None) -> SurfaceChamfer` — Chamfer between two intersecting surfaces: the ruled surface between the curves at the distances distance1 on s1 and distance2 on s2 (default distance1) from their intersection curve, measured within the surfaces (geodesics normal to the curve). Sides as for fillet_surfaces: +1 puts the chamfer on the side of the normal S_u x S_v. Returns a SurfaceChamfer (surface, edge, contact curves, p-curves on s1 / s2).

**Developable surfaces**

- `developability(surface: NurbsSurface, angular_tolerance: float = 1e-06) -> Developability` — Check whether a surface is developable: ruled along a parameter direction with the same normal all along every ruling. Reports the ruling direction, the largest turn of the normal along a ruling and the largest |Gaussian curvature|.
- `unroll(surface: NurbsSurface, tolerance: float = 1e-07, angular_tolerance: float = 1e-06) -> UnrolledSurface` — Unroll a developable surface (cylinder, cone, any ruled surface with constant normals along its rulings) into the plane without stretching. The pattern keeps the parameter domain: pattern(u, v) is the image of S(u, v), so curves on the surface map through their parameters. Raises GeometryError if the surface is not developable.

**Fitting**

- `interpolate_curve(points: array_like, degree: int = 3, method: Parametrization = Parametrization.Centripetal, end_derivatives: tuple[Vec3, Vec3] | None = None) -> NurbsCurve` — Global interpolation through points (n, 3); optional (C'(0), C'(1)).
- `approximate_curve(points: array_like, degree: int, n_ctrl: int, method: Parametrization = Parametrization.Centripetal) -> NurbsCurve` — Least-squares B-spline of the given `degree` with `n_ctrl` control points through the data `points` (shape (n, 3)); end points are interpolated, data parameters by `method`, knots by averaging (P&T Eq. 9.68).
- `approximate_curve_tol(points: array_like, degree: int, tolerance: float, method: Parametrization = Parametrization.Centripetal) -> NurbsCurve` — Least squares with the fewest control points meeting the tolerance.
- `interpolate_surface(points: array_like, degree_u: int = 3, degree_v: int = 3, method: Parametrization = Parametrization.Centripetal) -> NurbsSurface` — Interpolate a point grid (nu, nv, 3).
- `approximate_surface(points: array_like, degree_u: int, degree_v: int, ctrl_u: int, ctrl_v: int, method: Parametrization = Parametrization.Centripetal) -> NurbsSurface` — Least-squares surface with `ctrl_u` x `ctrl_v` control points to a point grid of shape (nu, nv, 3): first along u for every grid column, then along v.
- `fit_surface(points: array_like, num_u: int = 0, num_v: int = 0, degree_u: int = 3, degree_v: int = 3, smoothing: float = 0.001, iterations: int = 2, base: NurbsSurface | None = None, boundary: Sequence[NurbsCurve] = []) -> FittedSurface` — Fit a B-spline surface (domain [0, 1]^2) to unordered points (n, 3): least squares plus smoothing x thin-plate energy (relative weight; 1e-3 barely changes smooth data, 0.1-1 smooths noise), then `iterations` rounds of parameter correction. Parameters from the best-fit plane (height fields), from `base` (closest points), or from the Coons patch of four `boundary` curves [c0, c1, d0, d1], which are then interpolated exactly. num_u/num_v = 0: about sqrt(n / 4). Returns FittedSurface (surface, max_error, rms_error, parameters).
- `data_parameters(points: array_like, method: Parametrization = Parametrization.Centripetal) -> list[float]` — The parameters in [0, 1] assigned to the data points by the fitting functions for the given `method` (first 0, last 1).
- `fit_line(points: array_like) -> LineFit` — Best-fit line through points (n, 3) by orthogonal regression (minimal sum of squared perpendicular distances): centroid, unit direction, max and rms deviation.
- `fit_plane(points: array_like) -> PlaneFit` — Best-fit plane through points (n, 3) by orthogonal regression: plane through the centroid with unit normal, max and rms deviation.

**Analysis and properties**

- `(1) closest_point(curve: NurbsCurve, point: Vec3, tolerance: Tolerance = Tolerance()) -> CurvePoint`  
  `(2) closest_point(surface: NurbsSurface, point: Vec3, tolerance: Tolerance = Tolerance()) -> SurfacePoint` — (1) Global closest point on a curve (point projection / inversion). (2) Global closest point on a surface.
- `closest_point_local(curve: NurbsCurve, point: Vec3, u0: float, tolerance: Tolerance = Tolerance()) -> CurvePoint`  
  `closest_point_local(surface: NurbsSurface, point: Vec3, u0: float, v0: float, tolerance: Tolerance = Tolerance()) -> SurfacePoint` — Local Newton projection of `point` onto the curve (start parameter `u0`) or the surface (start `u0`, `v0`); fast, but converges to the nearest local minimum only. Use `closest_point` for the global minimum.
- `curvature(curve: NurbsCurve, u: float) -> float` — Curvature `|C' x C''| / |C'|^3`.
- `torsion(curve: NurbsCurve, u: float) -> float` — Torsion \tau(u) of the curve at `u` (zero for planar curves).
- `length(curve: NurbsCurve, rel_tol: float = 1e-10) -> float` — Total arc length of a curve.
- `arc_length(curve: NurbsCurve, u0: float, u1: float, rel_tol: float = 1e-10) -> float` — Arc length of the curve between the parameters `u0` and `u1` (adaptive Gauss–Kronrod quadrature to `rel_tol`).
- `parameter_at_length(curve: NurbsCurve, s: float, rel_tol: float = 1e-10) -> float` — Parameter u at which the arc length measured from the start of the domain equals `s` (Newton iteration on the arc length).
- `surface_curvature(surface: NurbsSurface, u: float, v: float) -> SurfaceCurvature` — Gaussian, mean and principal curvatures and directions.
- `area(surface: NurbsSurface, rel_tol: float = 1e-10) -> float`  
  `area(surface: NurbsSurface, u0: float, u1: float, v0: float, v1: float, rel_tol: float = 1e-10) -> float` — Area of an (untrimmed) surface, or of the parameter rectangle [`u0`, `u1`] x [`v0`, `v1`], by adaptive quadrature of `|S_u x S_v|`.
- `integrate(f: Callable[[float], float], a: float, b: float, rel_tol: float = 1e-10, breaks: Sequence[float] = [], abs_tol: float = 0.0) -> float` — Globally adaptive Gauss-Kronrod integral of f over [a, b].

**Curve analysis**

- `curve_curvature(curve: NurbsCurve, u: float) -> CurveCurvature` — Curvature, radius of curvature, unit tangent, principal normal and centre of the osculating circle at u.
- `curvature_graph(curve: NurbsCurve, scale: float = 1.0, samples_per_span: int = 16) -> CurvatureGraph` — Curvature graph: hairs from the curve points to point - scale * curvature vector (on the convex side, length proportional to the curvature), samples_per_span + 1 per knot span. Knots where the curvature jumps appear twice, with both values.
- `curve_direction(curve: NurbsCurve, plane_normal: Vec3 = (0, 0, 0), tolerance: float = 1e-06) -> CurveDirection` — Direction of a curve: end points, end tangents, and for a closed curve the sense of rotation seen from the tip of plane_normal (zero: the normal of the curve's own plane). curve.reversed() flips the direction.
- `is_closed(curve: NurbsCurve, tolerance: float = 1e-09) -> bool` — True if the end points of the curve agree within the tolerance.
- `is_planar(curve: NurbsCurve, tolerance: float = 1e-09) -> bool` — True if the curve lies in a plane within the tolerance.
- `(1) continuity(a: NurbsCurve, b: NurbsCurve, position: float = 1e-06, angle: float = 1e-06, curvature: float = 0.0001) -> CurveContinuity`  
  `(2) continuity(a: NurbsCurve, end_a: CurveEnd, b: NurbsCurve, end_b: CurveEnd, position: float = 1e-06, angle: float = 1e-06, curvature: float = 0.0001) -> CurveContinuity`  
  `(3) continuity(a: NurbsSurface, side: SurfaceSide, b: NurbsSurface, position: float = 1e-06, angle: float = 1e-06, curvature: float = 0.0001, samples: int = 33) -> SurfaceContinuity` — (1) Continuity of two curves at their nearest pair of ends: G0 / G1 / G2 (and C0 / C1 / C2) with the gap, the angle between the tangents and the curvatures on both sides. Tolerances: position (gap), angle (radians), curvature (relative). (2) Continuity at the joint of the given ends of a and b. (3) Continuity of surface a along its side with surface b: G0 / G1 / G2 with the largest gap, angle between the tangent planes and difference of the curvature tensors.
- `continuity_breaks(curve: NurbsCurve, position: float = 1e-06, angle: float = 1e-06, curvature: float = 0.0001) -> list[CurveBreak]` — Knots of a curve where it is not curvature continuous (corners, tangent joints), with the continuity there; for a closed curve also the seam.
- `find_duplicate_curves(curves: Sequence[NurbsCurve], tolerance: float = 1e-06) -> list[list[int]]` — Groups of geometrically identical curves (same shape within the tolerance, whatever the direction or parametrisation): lists of indices, the first being the one to keep.

**Surface analysis**

- `curvature_map(surface: object, samples_per_span: int = 8) -> CurvatureMap` — Curvature analysis of a surface or trimmed face: Gaussian, mean and principal curvatures on a grid with samples_per_span + 1 parameters per knot span, their ranges with the parameters where they are attained, and the smallest radius of curvature.
- `draft_angle(surface: NurbsSurface, u: float, v: float, pull_direction: Vec3) -> float` — Draft angle (radians) at (u, v): asin(N . d) for the unit normal N and the pull direction d; positive where the face can be withdrawn, negative for an undercut.
- `draft_angle_map(surface: object, pull_direction: Vec3, samples_per_span: int = 8) -> DraftAngleMap` — Draft angle analysis of a surface or trimmed face for a pull direction: angles (radians) on a grid and their range.
- `surface_direction(surface: NurbsSurface, tolerance: float = 1e-09) -> SurfaceDirection` — Directions of increasing u and v and the normal S_u x S_v at the middle of the domain, closed directions, poles and planarity. surface.transposed() flips the normal.
- `edge_continuity(shell: Shell, position: float = 1e-06, angle: float = 1e-06, curvature: float = 0.0001, samples: int = 9) -> list[SurfaceContinuity]` — Continuity across every edge of a shell (one SurfaceContinuity per edge): geometric >= 1 marks smooth (tangent) edges, -1 free and non-manifold edges.
- `find_duplicate_surfaces(surfaces: Sequence[NurbsSurface], tolerance: float = 1e-06) -> list[list[int]]` — Groups of geometrically identical surfaces (same point set within the tolerance, whatever the parametrisation or orientation): lists of indices, the first to keep.

**Intersections and sections**

- `(1) intersect(a: NurbsCurve, b: NurbsCurve, tolerance: Tolerance = Tolerance()) -> list[CurveCurveIntersection]`  
  `(2) intersect(curve: NurbsCurve, surface: NurbsSurface, tolerance: Tolerance = Tolerance()) -> list[CurveSurfaceIntersection]`  
  `(3) intersect(a: NurbsSurface, b: NurbsSurface, options: SurfaceIntersectionOptions = SurfaceIntersectionOptions()) -> SurfaceIntersection`  
  `(4) intersect(curve: NurbsCurve, plane: Plane, tolerance: Tolerance = Tolerance()) -> list[CurveSurfaceIntersection]` — (1) All intersection points of two curves, sorted by the parameter on a. (2) All intersection points of a curve and a surface, sorted by the curve parameter. (3) Surface/surface intersection: curves (3D + p-curves on both surfaces) and touching points. (4) Intersection points of a curve with a plane.
- `planar_section(surface: NurbsSurface, plane: Plane, options: SurfaceIntersectionOptions = SurfaceIntersectionOptions()) -> SurfaceIntersection` — Section of a surface with a plane (pcurve2 in the plane's in-plane coordinates).
- `planar_sections(surface: NurbsSurface, plane: Plane, spacing: float, count: int, options: SurfaceIntersectionOptions = SurfaceIntersectionOptions()) -> list[SurfaceIntersection]` — `count` planar sections of the surface with planes parallel to `plane`, at distances `k*spacing` along its normal (k = 0, ..., count-1).

**Trimming and splitting**

- `rectangle_loop(u0: float, u1: float, v0: float, v1: float) -> TrimLoop` — Counter-clockwise rectangular trim loop [`u0`, `u1`] x [`v0`, `v1`] in parameter space.
- `circle_loop(cu: float, cv: float, r: float) -> TrimLoop` — Counter-clockwise exact circular trim loop with centre (`cu`, `cv`) and radius `r` in parameter space.
- `polygon_loop(points: Sequence[Sequence[float] (2)]) -> TrimLoop` — Closed polygonal trim loop through the (u, v) points (closed automatically). Use counter-clockwise order for an outer loop and clockwise order for a hole.
- `(1) split(surface: NurbsSurface, cutter: NurbsSurface, tolerance: Tolerance = Tolerance()) -> list[TrimmedSurface]`  
  `(2) split(surface: NurbsSurface, pcurves: Sequence[NurbsCurve], tolerance: Tolerance = Tolerance()) -> list[TrimmedSurface]`  
  `(3) split(face: TrimmedSurface, cuts: Sequence[NurbsCurve], tolerance: Tolerance = Tolerance()) -> list[TrimmedSurface]`  
  `(4) split(shell: Shell, plane: Plane, tolerance: float = 1e-06, cut_label: str = 'cut') -> tuple[Shell, Shell]` — (1) Split a surface along its intersection curves with another surface into trimmed regions. (2) Split a surface into trimmed regions along parameter-space curves (u, v, 0). Closed curves become holes and new regions; open curves run from boundary to boundary and may cross. (3) Split a trimmed face along parameter-space curves (general planar arrangement): cuts are clipped to the face, may end on the trim loops or on other cuts and may cross each other. (4) Split a closed shell by a plane -> (below, above): two closed shells, below opposite to the plane normal. The new faces in the plane get `cut_label` (e.g. for a symmetry boundary condition); a part is an empty shell if the plane misses the solid.

**Booleans and shells**

- `boolean(a: Shell, b: Shell, op: BooleanOp, tolerance: float = 1e-06) -> Shell` — Boolean union / intersection / difference of two closed shells.
- `point_in_shell(shell: Shell, point: Vec3, tolerance: float = 1e-06) -> bool` — True if the point lies inside the closed shell (ray casting).

**Solid operations**

- `cap(shell: Shell, tolerance: float = 1e-06, label: str = 'cap') -> Shell` — Close the planar openings of a shell with planar faces trimmed exactly by the opening's edges (coplanar nested openings: faces with holes). The new faces get `label`; the result is sewn and oriented outwards. Raises GeometryError for non-planar openings.
- `thicken(shell: Shell, thickness: float, symmetric: bool = False, tolerance: float = 1e-06, wall_label: str = 'wall') -> Shell` — Give an open shell a wall thickness: faces offset along S_u x S_v (negative: against it; symmetric: +-thickness/2), free edges closed by ruled side walls labelled `wall_label`, offset faces labelled '<label>_offset'. The faces must join smoothly; raises GeometryError at creases or where the thickness reaches a radius of curvature.
- `hollow(shell: Shell, thickness: float, open_labels: Sequence[str] = [], tolerance: float = 1e-06) -> Shell` — Hollow a smooth closed shell to a wall of the given thickness: an inner void, or - with open_labels - removing the faces with these labels and thickening the rest inwards.
- `(1) make_pipe(path: NurbsCurve, radius: float, end_radius: float = -1.0, tolerance: float = 1e-06) -> Shell`  
  `(2) make_pipe(path: NurbsCurve, radius: Callable[[float], float], tolerance: float = 1e-06) -> Shell` — (1) Solid pipe along a path: circular cross-section (radius varying linearly to `end_radius` if given), planar end faces. Faces labelled 'wall', 'inlet', 'outlet'; a closed path gives a ring. (2) Solid pipe whose radius is the function radius(t) of the normalised path parameter t in [0, 1] (smooth and positive).
- `remove_holes(shell: Shell, max_size: float, tolerance: float = 1e-06) -> Shell` — Defeaturing: fill the holes (through or blind) of a closed shell whose loops are smaller than max_size across; the faces of the removed features are deleted.
- `merge_faces(shell: Shell, tolerance: float = 1e-06) -> Shell` — Defeaturing: merge adjacent faces with the same label lying on the same surface (pieces of one split face, or coplanar planar faces with the same normal) into one face.
- `edges_between(shell: Shell, label_a: str, label_b: str) -> list[int]` — Indices of the edges shared by a face labelled `label_a` and one labelled `label_b`.
- `fillet_edges(shell: Shell, edges: Sequence[int], radius: float, label: str = 'fillet', tolerance: float = 0.0, sizes: Sequence[tuple[Vec3, float]] = []) -> Shell` — Round edges of a closed shell with a rolling-ball fillet (convex edges lose, concave edges gain material). Edges between the same two faces are filleted as one chain. Where chains meet at a vertex of three faces: tangent chains continue each other, two chains are mitred (trimmed by each other), three chains get a spherical corner patch; a chain ending alone is trimmed by the third face. Variable radius: `sizes` = [(point, radius), ...] at points on the selected edges, smooth (C2) in between, constant beyond; chains without points use `radius`. New faces are labelled `label`.
- `chamfer_edges(shell: Shell, edges: Sequence[int], distance: float, label: str = 'chamfer', tolerance: float = 0.0, sizes: Sequence[tuple[Vec3, float]] = []) -> Shell` — Chamfer edges of a closed shell: a ruled cut between the curves at (geodesic) `distance` from the edge on both faces. Grouping, vertices, variable `sizes` (point, distance) and restrictions as for fillet_edges; a corner of three chamfered edges (between planar faces) gets a planar triangle.

**Slicing and hatching**

- `section(shell: Shell, plane: Plane, options: SectionOptions = SectionOptions()) -> PlanarSection` — Section of a closed shell with a plane: closed loops oriented with the solid on the left (outer counter-clockwise, holes clockwise in the plane coordinates).
- `critical_heights(shell: Shell, direction: Vec3, tolerance: Tolerance = Tolerance()) -> list[float]` — Sorted heights along `direction` at which the topology of the planar sections can change.
- `is_extrusion_band(shell: Shell, direction: Vec3, h0: float, h1: float, tolerance: Tolerance = Tolerance()) -> bool` — True if all faces meeting the band h0 < p . direction < h1 are parallel to `direction` there.
- `hatch(loops: object, direction: Sequence, spacing: float, *, phase: float = 0.0, contours: Sequence[float] = [], min_length: float = 0.0, tolerance: float = 0.0, contour_tolerance: float = 0.0) -> tuple[list[HatchSegment], list[HatchContour]]` — Hatch and contour lines of the region bounded by closed planar loops (region on the left). `loops`: a PlanarSection, its loops, or closed chains of curves in the xy plane (a single closed curve counts as a loop). Lines n . p = phase + k * spacing along `direction` (2 or 3 components, z ignored), n = (-d.y, d.x). Returns (segments, contours): HatchSegment list sorted by line and position, and one HatchContour per loop and inward offset distance in `contours`.

**Silhouettes and drawings**

- `silhouette(object: object, direction: Vec3 | None = None, eye: Vec3 | None = None, tolerance: float = 1e-06, resolution: int = 8) -> list[SilhouetteCurve]` — Silhouette curves of a surface, trimmed face or shell: where the normal is perpendicular to the view ray. Parallel view along `direction` (from the viewer into the scene) or perspective view from the point `eye`. Curves are fitted to `tolerance`; `resolution` is the number of search cells per knot span. Faces seen edge-on have none.
- `visibility(curves: Sequence[NurbsCurve], occluders: Sequence[Shell], direction: Vec3 | None = None, eye: Vec3 | None = None, tolerance: float = 1e-06, samples: int = 64) -> list[CurvePiece]` — Hidden-line removal: split every curve where its visibility against the shells changes and classify the pieces (CurvePiece.visible). Parallel view along `direction` or perspective view from `eye`. Curves lying on a face are not hidden by that face.
- `make2d(shells: Sequence[Shell], direction: Vec3, s_axis: Vec3 = (0, 0, 0), origin: Vec3 = (0, 0, 0), tolerance: float = 1e-06, samples: int = 64, resolution: int = 8, join_tolerance: float = 0.0, hidden: bool = True, tangent_edges: bool = False) -> Drawing` — 2D drawing of shells for a parallel view along `direction` (from the viewer into the scene): all edges and silhouettes, split by visibility and projected onto the plane through `origin` normal to the direction, in the coordinates (s, t) with s along `s_axis` (zero: automatic) and t = normal x s. Returns a Drawing with the visible and the hidden lines; duplicates are removed, visible lines taking precedence. tangent_edges: also draw smooth edges (seams, borders of fillets).
- `outline(shells: Sequence[Shell], direction: Vec3, s_axis: Vec3 = (0, 0, 0), origin: Vec3 = (0, 0, 0), tolerance: float = 1e-06, resolution: int = 8, join_tolerance: float = 0.0) -> Outline` — Contour of shells seen along `direction`: the boundary of the union of their parallel projections, as closed loops in the plane coordinates of make2d (Outline.loops, a region for region_boolean / region_area: outer contours counter-clockwise, holes clockwise).

**Optical surfaces**

- `asphere_surface(asphere: Asphere, r_max: float, tolerance: float = 1e-07) -> NurbsSurface` — Surface of revolution; exact for pure conics.
- `asphere_profile(asphere: Asphere, r_max: float, tolerance: float = 1e-07) -> NurbsCurve` — The meridian profile z(r), r in [0, `r_max`], of the asphere as a planar curve in the xz plane: exact conic (rational) without polynomial terms, otherwise a cubic spline within `tolerance`. `asphere_surface` revolves it.
- `zernike(j: int, indexing: ZernikeIndexing, rho: float, theta: float) -> float` — Value of the Zernike polynomial Z_j(\rho, \theta) in the given indexing (j ≥ 1).
- `zernike_nm(j: int, indexing: ZernikeIndexing) -> tuple[int, int]` — Radial degree n and azimuthal frequency m (signed: m < 0 for sine terms) of the Zernike term j in the given indexing.
- `zernike_surface(zernike: ZernikeSag, half_width: float = 0.0, tolerance: float = 1e-07) -> NurbsSurface` — The Zernike sag (base asphere plus Zernike terms) as a graph surface over [-h, h]^2 with u = x, v = y; `half_width` h = 0 uses `norm_radius`. Fitted until the sag error is at most `tolerance` at the check points.
- `xy_polynomial_surface(poly: XYPolynomialSag, half_width: float = 0.0, tolerance: float = 1e-07) -> NurbsSurface` — The XY-polynomial sag as a graph surface over [-h, h]^2 (u = x, v = y); `half_width` h = 0 uses `norm_radius`. Fitted to `tolerance`.
- `sag_surface(f: Callable[[float, float], float], x0: float, x1: float, y0: float, y1: float, degree: int = 3, tolerance: float = 1e-07, max_spans: int = 256) -> NurbsSurface` — Graph surface z = f(x, y) over a rectangle, fitted to the tolerance.

**STEP import and export**

- `write_step(path: object, shells: object, unit: str = 'mm', name: str = 'nurbforge', author: str = '', organization: str = '', timestamp: str = '', uncertainty: float = 1e-06, names: Sequence[str] = []) -> None` — Write a Shell or a list of Shells to a STEP file (AP214, native writer). Closed shells become solids, open shells surface models. `unit` ('mm', 'cm', 'm', 'inch', 'um') declares what one model unit means; coordinates are written unchanged. `names`: one name per shell (solid names, read back by read_step_named).
- `to_step(shells: object, unit: str = 'mm', name: str = 'nurbforge', author: str = '', organization: str = '', timestamp: str = '', uncertainty: float = 1e-06, names: Sequence[str] = []) -> str` — The STEP file contents (see write_step) as a string.
- `read_step(path: object, *, unit: str = 'mm', tolerance: float = 1e-06, skip_unsupported: bool = False) -> list[Shell]` — Read the B-rep solids and surface models of a STEP file (native reader) as a list of Shells, converted to `unit` ('mm', 'cm', 'm', 'inch', 'um'). Analytic surfaces (plane, cylinder, cone, sphere, torus) become exact NURBS patches; trim loops are rebuilt in parameter space. `tolerance` is the sewing tolerance. Unsupported entities raise GeometryError unless `skip_unsupported` (the faces are then skipped with a warning).
- `from_step(text: str, *, unit: str = 'mm', tolerance: float = 1e-06, skip_unsupported: bool = False) -> list[Shell]` — Read shells from STEP file contents (see read_step).
- `read_step_named(path: object, *, unit: str = 'mm', tolerance: float = 1e-06, skip_unsupported: bool = False) -> list` — As read_step, as a list of (name, Shell): the solid name (write_step(..., names=...)), else the product name.
- `from_step_named(text: str, *, unit: str = 'mm', tolerance: float = 1e-06, skip_unsupported: bool = False) -> list` — As from_step, as a list of (name, Shell).

**Transforms and patterns**

- `linear_pattern(step: Vec3, count: int) -> list[Transform]` — `count` translations by 0, `step`, 2 `step`, ... (the first is the identity).
- `grid_pattern(step_a: Vec3, count_a: int, step_b: Vec3, count_b: int) -> list[Transform]` — `count_a` x `count_b` translations `i*step_a + j*step_b`.
- `polar_pattern(axis_point: Vec3, axis_dir: Vec3, count: int, total_angle: float) -> list[Transform]` — `count` rotations about the axis through `axis_point` along `axis_dir`, evenly spread over `total_angle` radians (a full turn: angles k 2π/count; otherwise the last copy lands at `total_angle`).
- `apply_pattern(geometry: G, transforms: Iterable[Transform]) -> list[G]` — Return one transformed copy of `geometry` per transform.

**Knot vectors**

- `clamped_uniform_knots(degree: int, n_ctrl: int) -> list[float]` — Clamped knot vector on [0, 1] with uniform interior knots for `n_ctrl` control points of the given degree (length `n_ctrl + degree + 1`).
- `unique_knots(knots: Sequence[float]) -> tuple[list[float], list[int]]` — Distinct knot values and multiplicities (gmsh/OCC format).
- `validate_knots(knots: Sequence[float], degree: int, n_ctrl: int) -> None` — Check a knot vector for `n_ctrl` control points of the given degree (length, monotonicity, multiplicities); raises `GeometryError` with the reason, returns `None` otherwise.

**Vector algebra**

- `dot(arg0: Vec3, arg1: Vec3) -> float` — Dot product of two vectors.
- `cross(arg0: Vec3, arg1: Vec3) -> Vec3` — Cross product of two vectors.
- `norm(arg0: Vec3) -> float` — Euclidean length of a vector.
- `distance(arg0: Vec3, arg1: Vec3) -> float` — Euclidean distance between two points.

**Help**

- `reference() -> str` — The compact API reference for language models and coding agents (Markdown): conventions, a worked example and every public function, class and LLM tool with signature and description. `python -m nurbforge` prints it.

## Classes of `nurbforge`

### `Asphere`

Rotationally symmetric sag: conic + even (r^4, r^6, ...) + odd (r, r^2, ...) terms. Rotationally symmetric aspheric sag z(r) = c r^2 / (1 + \sqrt{1 - (1+k) c^2 r^2}) + \sum_i A_{2i+2} r^{2i+2} + \sum_i B_i r^i with curvature c, conic constant k, even coefficients A4, A6, ... and odd coefficients B1, B2, .... Constructor arguments are keywords: `curvature`, `conic`, `even`, `odd`.

- `Asphere(curvature: float = 0.0, conic: float = 0.0, even: Sequence[float] = [], odd: Sequence[float] = [])` — Create an asphere; without `even`/`odd` terms it is an exact conic.
- `.conic: float` (settable) — Conic constant k (0 sphere, -1 paraboloid, < -1 hyperboloid, (-1, 0) prolate and > 0 oblate ellipsoid).
- `.curvature: float` (settable) — Vertex curvature c = 1/R (R > 0: centre of curvature on +z).
- `.even: list[float]` (settable) — Even polynomial coefficients A4, A6, A8, ... of r^4, r^6, r^8, …
- `.is_conic: bool` — True if there are no polynomial terms (the surface is an exact conic).
- `.max_radius() -> float` — Largest r for which the conic term is defined (infinite for k ≤ -1 or c = 0).
- `.odd: list[float]` (settable) — Odd/general polynomial coefficients B1, B2, B3, ... of r, r^2, r^3, …
- `.sag(r: float) -> float` — Sag z at radius r.

### `CurvatureGraph`

Curvature graph (comb) of a curve.

- `.curvatures: list[float]` — Curvature \kappa at the samples.
- `.max_curvature: float` — Largest curvature among the samples.
- `.max_parameter: float` — Parameter of the sample with the largest curvature.
- `.parameters: list[float]` — Parameters of the samples (a knot with a curvature jump appears twice).
- `.points: list[Vec3]` — curve points (feet of the hairs)
- `.tips: list[Vec3]` — ends of the hairs; in order, the graph line

### `CurvatureMap`

Curvatures of a surface on a grid; arrays of shape (len(u), len(v)).

- `.gaussian: ndarray` — Gaussian curvature K = k_1 k_2, array of shape (len(u), len(v)).
- `.gaussian_range: SurfaceMapRange` — `SurfaceMapRange` of the Gaussian curvature.
- `.inside: numpy.NDArray[numpy.bool]` — Boolean array: samples inside the trimmed region (all true for a surface).
- `.k_max: ndarray` — Larger principal curvature.
- `.k_max_range: SurfaceMapRange` — `SurfaceMapRange` of the larger principal curvature.
- `.k_min: ndarray` — Smaller principal curvature.
- `.k_min_range: SurfaceMapRange` — `SurfaceMapRange` of the smaller principal curvature.
- `.mean: ndarray` — Mean curvature H = (k_1 + k_2)/2, sign relative to S_u × S_v.
- `.mean_range: SurfaceMapRange` — `SurfaceMapRange` of the mean curvature.
- `.min_radius: float` — Smallest radius of curvature among the samples (inf for a plane).
- `.u: list[float]` — Parameters u of the grid.
- `.v: list[float]` — Parameters v of the grid.

### `CurveBreak`

Parameter of a curve with its continuity there.

- `.continuity: CurveContinuity` — `CurveContinuity` at the parameter, from the one-sided derivatives.
- `.parameter: float` — Parameter of the knot (domain end for the seam of a closed curve).

### `CurveChamfer`

Chamfer (straight cut) between two curves.

- `.line: NurbsCurve` — straight segment from point_a to point_b
- `.param_a: float` — parameter of the chamfer point on a
- `.param_b: float` — parameter of the chamfer point on b
- `.point_a: Vec3` — chamfer point on a (line start)
- `.point_b: Vec3` — chamfer point on b (line end)
- `.trimmed_a: NurbsCurve | None` — kept part of a ending at point_a (None if trim=False)
- `.trimmed_b: NurbsCurve | None` — kept part of b starting at point_b (None if trim=False)

### `CurveContinuity`

Continuity of curves at a joint.

- `.curvature_a: float` — Curvature on the side of the first curve.
- `.curvature_b: float` — Curvature on the side of the second curve.
- `.curvature_difference: float` — Norm of the difference of the two curvature vectors.
- `.gap: float` — Distance of the two end points.
- `.geometric: int` — -1 no contact, 0 position (G0), 1 tangent (G1), 2 curvature (G2)
- `.parametric: int` — -1, 0 (C0), 1 (C1), 2 (C2)
- `.tangent_angle: float` — radians

### `CurveCurvature`

Curvature of a curve at a parameter, with the osculating circle.

- `.center: Vec3` — centre of the osculating circle
- `.curvature: float` — Curvature \kappa ≥ 0.
- `.normal: Vec3` — unit principal normal, towards the centre of curvature
- `.point: Vec3` — The curve point C(u).
- `.radius: float` — 1 / curvature (inf where it vanishes)
- `.tangent: Vec3` — unit tangent
- `.vector: Vec3` — curvature vector: curvature * normal

### `CurveCurveIntersection`

One intersection of two curves: parameters `s` on the first and `t` on the second curve, the point and the remaining distance between the curves there.

- `.distance: float` — Distance between the two curve points (0 up to the tolerance).
- `.point: Vec3` — Intersection point (midpoint of the two curve points).
- `.s: float` — Parameter on the first curve.
- `.t: float` — Parameter on the second curve.

### `CurveDirection`

Direction of a curve and sense of rotation of a closed curve.

- `.closed: bool` — True if the end points agree within the tolerance.
- `.end: Vec3` — End point.
- `.end_tangent: Vec3` — Unit tangent at the end, in the direction of travel.
- `.orientation: int` — +1 counter-clockwise, -1 clockwise seen from the tip of plane_normal; 0 open
- `.planar: bool` — True if the curve lies in a plane within the tolerance.
- `.plane_normal: Vec3` — Normal the sense of rotation refers to: the given one, or that of the curve's plane with its largest component positive.
- `.signed_area: float` — Enclosed area of the projection along the plane normal, positive for counter-clockwise curves.
- `.start: Vec3` — Start point.
- `.start_tangent: Vec3` — Unit tangent at the start.

### `CurveFillet`

Circular fillet between two planar curves.

- `.arc: NurbsCurve` — exact rational arc from point_a to point_b
- `.center: Vec3` — Centre of the fillet arc.
- `.param_a: float` — parameter of the tangent point on a
- `.param_b: float` — parameter of the tangent point on b
- `.point_a: Vec3` — Tangent point on curve a, `a(param_a)` = start of the arc.
- `.point_b: Vec3` — Tangent point on curve b, `b(param_b)` = end of the arc.
- `.trimmed_a: NurbsCurve | None` — kept part of a ending at point_a (None if trim=False)
- `.trimmed_b: NurbsCurve | None` — kept part of b starting at point_b (None if trim=False)

### `CurvePiece`

Piece of a curve with constant visibility.

- `.curve: NurbsCurve` — The part of the source curve on [t0, t1] (same parameter).
- `.source: int` — index of the source curve
- `.t0: float` — Start of the parameter interval on the source curve.
- `.t1: float` — End of the parameter interval on the source curve.
- `.visible: bool` — True if the piece is visible.

### `CurvePoint`

Result of a closest-point query on a curve.

- `.distance: float` — Distance from the query point.
- `.point: Vec3` — The closest point on the curve.
- `.u: float` — Parameter of the closest point.

### `CurveSurfaceIntersection`

One intersection of a curve with a surface or plane: curve parameter `t`, surface parameters `u`, `v`, point and remaining distance.

- `.distance: float` — Distance between curve and surface point (0 up to the tolerance).
- `.point: Vec3` — Intersection point.
- `.t: float` — Parameter on the curve.
- `.u: float` — u parameter on the surface (in-plane coordinate s for a plane).
- `.v: float` — v parameter on the surface (in-plane coordinate t for a plane).

### `Developability`

Result of developability().

- `.developable: bool` — True if the surface is ruled and the normal is the same along every ruling (within the angular tolerance).
- `.max_gaussian_curvature: float` — Largest |K| at the sample points (0 for a developable surface).
- `.max_twist: float` — largest angle (rad) between normals along one ruling
- `.rulings: RulingDirection` — `RulingDirection` of the straight parameter lines.

### `DraftAngleMap`

Draft angles of a surface on a grid; arrays of shape (len(u), len(v)).

- `.angle: ndarray` — Draft angles in radians, array of shape (len(u), len(v)).
- `.inside: numpy.NDArray[numpy.bool]` — Boolean array: samples inside the trimmed region.
- `.range: SurfaceMapRange` — `SurfaceMapRange` of the angle over the samples inside.
- `.u: list[float]` — Parameters u of the grid.
- `.v: list[float]` — Parameters v of the grid.

### `Drawing`

Visible and hidden lines of shells in a drawing plane. Result of `make2d`: the visible and the hidden lines in a plane with the frame (origin, s_axis, t_axis, normal); a point is origin + s s_axis + t t_axis.

- `.hidden: list[DrawingCurve]` — Hidden lines (`DrawingCurve`); empty with `hidden=False`.
- `.normal: Vec3` — Unit normal of the drawing plane, pointing to the viewer.
- `.origin: Vec3` — Origin of the plane coordinates.
- `.s_axis: Vec3` — First in-plane unit axis s (to the right, seen from the viewer).
- `.t_axis: Vec3` — Second in-plane unit axis t = normal x s (upwards).
- `.visible: list[DrawingCurve]` — Visible lines (`DrawingCurve`).

### `DrawingCurve`

One line of a drawing.

- `.curve: NurbsCurve` — in the drawing plane: (s, t, 0)
- `.curve3d: NurbsCurve` — the same piece in 3D
- `.index: int` — edge index, or face index of a silhouette
- `.kind: DrawingCurveKind` — `DrawingCurveKind`: edge of a shell or silhouette of a face.
- `.shell: int` — index of the source shell
- `.visible: bool` — True for visible lines.

### `FittedSurface`

Result of fit_surface.

- `.max_error: float` — largest distance of a point from the surface
- `.parameters: ndarray` — (n, 2) array: (u, v) of each point's foot point
- `.rms_error: float` — rms of the point distances
- `.surface: NurbsSurface` — the fitted surface

### `GeometryError` (exception, subclass of ValueError)

Raised for invalid input or degenerate geometry (invalid knot vectors, zero-length axes, tolerances that cannot be met, open shells where closed ones are required, ...). The message names the cause. Subclass of `ValueError`.

### `HatchContour`

Closed contour: boundary of the region shrunk by `distance`.

- `.curves: list[NurbsCurve]` — Consecutive pieces of the closed contour (offset curves and arcs at reflex corners), in the xy plane.
- `.distance: float` — Inset distance of this contour from the boundary.
- `.index: int` — Index of the distance in `contours`.
- `.loop: int` — First of `loops`.
- `.loops: list[int]` — Loops it was offset from (merged contours: several).

### `HatchSegment`

Hatch segment inside the region bounded by the loops.

- `.end: Vec3` — End point (x, y, 0).
- `.end_loop: int` — Index of the loop on which the segment ends.
- `.length: float` — Length of the segment.
- `.line: int` — Index k of the line n . p = phase + k * spacing.
- `.region: int` — Outer loop of the region containing the segment.
- `.start: Vec3` — Start point (x, y, 0); start lies before end along the hatch direction.
- `.start_loop: int` — Index of the loop on which the segment starts.

### `IntersectionCurve`

One branch of a surface/surface intersection: the traced points, the fitted 3D curve and its p-curves (u, v, 0) on both surfaces, which share the curve's parameter.

- `.closed: bool` — True if the branch is a closed loop.
- `.curve: NurbsCurve` — Fitted 3D intersection curve.
- `.pcurve1: NurbsCurve` — Curve (u1, v1, 0) in the parameter domain of the first surface.
- `.pcurve2: NurbsCurve` — Curve (u2, v2, 0) in the parameter domain of the second surface (in-plane coordinates for `planar_section`).
- `.points: list[SurfaceSurfacePoint]` — The traced points with their parameters on both surfaces.

### `LineFit`

Best-fit line (orthogonal regression).

- `.direction: Vec3` — unit direction
- `.max_deviation: float` — largest perpendicular distance of a point from the line
- `.point: Vec3` — centroid of the points
- `.rms: float` — root mean square of the perpendicular distances

### `MassProperties`

Volume, centroid and inertia of a solid (unit density).

- `.centroid: Vec3` — Centre of mass (unit density).
- `.inertia: ndarray` — Inertia tensor about the centroid (3 x 3): diagonal moments, off-diagonal -products.
- `.principal_axes: ndarray` — Unit axes as rows, in the order of principal_moments.
- `.principal_moments: tuple[float, float, float]` — Principal moments of inertia about the centroid (eigenvalues of `inertia`, ascending).
- `.volume: float` — Enclosed volume.

### `NurbsCurve`

Non-uniform rational B-spline curve: degree p, knot vector of length n + p + 2, control points P_i and weights w_i > 0. Immutable; every modifying operation returns a new curve. Evaluate with `c(u)`.

- `NurbsCurve(degree: int, knots: Sequence[float], control_points: array_like, weights: Sequence[float] = [])` — Create a curve from its degree, the full knot vector, the control points (array of shape (n, 3)) and optional weights (empty: non-rational). Validated; raises `GeometryError`.
- `.clamped() -> NurbsCurve` — Same curve with a clamped knot vector.
- `.control_points: ndarray` — Control points as an array of shape (n, 3).
- `.decompose_bezier() -> list[NurbsCurve]` — Bezier segments, each on its original parameter interval.
- `.degree: int` — Polynomial degree p.
- `.derivatives(u: float, order: int = 1) -> list[Vec3]` — The point and its derivatives up to `order` at u: `[C(u), C'(u), C''(u), ...]` (rational derivatives).
- `.domain: tuple[float, float]` — Parameter range (u_min, u_max).
- `.elevate_degree(times: int = 1) -> NurbsCurve` — Degree elevation (shape-preserving); result is clamped.
- `.evaluate(u: float) -> Vec3` — Point C(u); same as `c(u)`.
- static `NurbsCurve.from_points(degree: int, control_points: array_like) -> NurbsCurve` — Clamped uniform B-spline with the given degree and control points (not an interpolation; see `interpolate_curve`).
- `.insert_knot(u: float, times: int = 1) -> NurbsCurve` — Insert the knot u `times` times (Boehm, P&T A5.1); the geometry is unchanged.
- `.is_clamped: bool` — True if the first and last knots have multiplicity p + 1 (the curve starts and ends at its end control points).
- `.is_rational: bool` — True if the weights are not all equal.
- `.knots: list[float]` — The full knot vector.
- `.reduce_degree(tolerance: float = 1e-10) -> NurbsCurve` — Degree reduction by one within tolerance (raises GeometryError otherwise).
- `.refine_knots(knots: Sequence[float]) -> NurbsCurve` — Insert all given knots at once (shape-preserving).
- `.remove_knot(u: float, times: int = 1, tolerance: float = 1e-10) -> NurbsCurve` — Remove the knot u up to `times` times as long as the curve moves by at most `tolerance`; raises `GeometryError` if not even one removal is possible (see `try_remove_knot`).
- `.reversed() -> NurbsCurve` — The same curve traversed in the opposite direction (same domain).
- `.sample(n: int = 100) -> ndarray` — `n` points at equally spaced parameters, as an array of shape (n, 3).
- `.split(u: float) -> tuple[NurbsCurve, NurbsCurve]` — Split at u -> (left, right); both keep their parameter intervals.
- `.sub_curve(u0: float, u1: float) -> NurbsCurve` — The piece between the parameters `u0` < `u1` (exact, by knot insertion), keeping their parameters.
- `.tangent(u: float) -> Vec3` — Unit tangent at u.
- `.transformed(transform: Transform) -> NurbsCurve` — The curve with every control point transformed (exact for affine maps).
- `.try_remove_knot(u: float, times: int = 1, tolerance: float = 1e-10) -> tuple[NurbsCurve, int]` — Remove knot u up to `times` times within tolerance -> (curve, removed).
- `.weights: list[float]` — The weights (all 1 for a non-rational curve).
- `__call__(u: float) -> Vec3` — Call operator.
- `__rmul__(arg0: Transform) -> NurbsCurve` — Operator `b * a`.

### `NurbsSurface`

Tensor-product NURBS surface of degrees (p, q) with a net of `num_u` x `num_v` control points, stored row-major (index i*num_v + j, i along u). Immutable. Evaluate with `s(u, v)`.

- `NurbsSurface(degree_u: int, degree_v: int, knots_u: Sequence[float], knots_v: Sequence[float], num_u: int, num_v: int, control_points: array_like, weights: Sequence[float] = [])` — Control points row-major: index = i * num_v + j (i along u).
- `.clamped() -> NurbsSurface` — The same surface with clamped knot vectors in u and v.
- `.control_points: ndarray` — Control net as an array of shape (num_u * num_v, 3), row-major (index i * num_v + j).
- `.decompose_bezier() -> list[NurbsSurface]` — Bezier patches in row-major order (u span outer), original parameter rectangles.
- `.degree_u: int` — Degree p in u.
- `.degree_v: int` — Degree q in v.
- `.derivatives(u: float, v: float, order: int = 1) -> list[list[Vec3]]` — Partial derivatives D[k][l] = d^(k+l)S/du^k dv^l for k + l <= order.
- `.domain_u: list[float] (2)` — Parameter range [u_min, u_max].
- `.domain_v: list[float] (2)` — Parameter range [v_min, v_max].
- `.elevate_degree_u(times: int = 1) -> NurbsSurface` — Raise the degree in u by `times` (geometry unchanged).
- `.elevate_degree_v(times: int = 1) -> NurbsSurface` — Raise the degree in v by `times` (geometry unchanged).
- `.evaluate(u: float, v: float) -> Vec3` — Point S(u, v); same as `s(u, v)`.
- `.first_derivatives(u: float, v: float) -> list[Vec3] (3)` — `[S, S_u, S_v]` at (u, v).
- `.iso_curve_u(u: float) -> NurbsCurve` — The curve v \mapsto S(u, v) at fixed u (exact).
- `.iso_curve_v(v: float) -> NurbsCurve` — The curve u \mapsto S(u, v) at fixed v (exact).
- `.knots_u: list[float]` — Knot vector in u.
- `.knots_v: list[float]` — Knot vector in v.
- `.normal(u: float, v: float) -> Vec3` — Unit normal along `S_u x S_v`; at degenerate points (poles of revolved surfaces) the limit normal.
- `.num_u: int` — Number of control points in u.
- `.num_v: int` — Number of control points in v.
- `.reduce_degree_u(tolerance: float = 1e-10) -> NurbsSurface` — Lower the degree in u by one if the surface moves by at most `tolerance`; raises `GeometryError` otherwise.
- `.reduce_degree_v(tolerance: float = 1e-10) -> NurbsSurface` — Lower the degree in v by one if the surface moves by at most `tolerance`; raises `GeometryError` otherwise.
- `.refine_knots_u(knots: Sequence[float]) -> NurbsSurface` — Insert all given knots in u at once (P&T A5.5); geometry unchanged.
- `.refine_knots_v(knots: Sequence[float]) -> NurbsSurface` — Insert all given knots in v at once; geometry unchanged.
- `.remove_knot_u(u: float, times: int = 1, tolerance: float = 1e-10) -> NurbsSurface` — Remove the u knot up to `times` times within `tolerance` (common count over all rows); raises if none can be removed.
- `.remove_knot_v(v: float, times: int = 1, tolerance: float = 1e-10) -> NurbsSurface` — Remove the v knot up to `times` times within `tolerance`; raises if none can be removed.
- `.sample(nu: int = 20, nv: int = 20) -> ndarray` — Grid of `nu` x `nv` points at equally spaced parameters, array of shape (nu, nv, 3).
- `.split_u(u: float) -> tuple[NurbsSurface, NurbsSurface]` — Split at u into two surfaces (exact).
- `.split_v(v: float) -> tuple[NurbsSurface, NurbsSurface]` — Split at v into two surfaces (exact).
- `.sub_surface(u0: float, u1: float, v0: float, v1: float) -> NurbsSurface` — The patch [`u0`, `u1`] x [`v0`, `v1`] (exact), keeping its parameters.
- `.transformed(transform: Transform) -> NurbsSurface` — The surface with every control point transformed.
- `.transposed() -> NurbsSurface` — The same surface with u and v exchanged (the normal flips).
- `.try_remove_knot_u(u: float, times: int = 1, tolerance: float = 1e-10) -> tuple[NurbsSurface, int]` — Remove the u knot up to `times` times within `tolerance`; returns (surface, number actually removed).
- `.try_remove_knot_v(v: float, times: int = 1, tolerance: float = 1e-10) -> tuple[NurbsSurface, int]` — Remove the v knot up to `times` times within `tolerance`; returns (surface, number actually removed).
- `.weights: list[float]` — Weights in net order (all 1 for a non-rational surface).
- `__call__(u: float, v: float) -> Vec3` — Call operator.
- `__rmul__(arg0: Transform) -> NurbsSurface` — Operator `b * a`.

### `Outline`

Contour of shells seen along a direction. Result of `outline`: closed loops in a plane with the frame (origin, s_axis, t_axis, normal).

- `.loops: list[list[NurbsCurve]]` — Region bounded by the contour: loops of curves (s, t, 0), outer contours counter-clockwise, holes clockwise.
- `.normal: Vec3` — Unit normal of the plane, pointing to the viewer.
- `.origin: Vec3` — Origin of the plane coordinates.
- `.s_axis: Vec3` — First in-plane unit axis s.
- `.t_axis: Vec3` — Second in-plane unit axis t = normal x s.

### `PlanarSection`

Result of `section`: closed loops in a plane with the frame (origin, s_axis, t_axis, normal). Curves are given in the plane coordinates (s, t, 0).

- `.loops: list[SectionLoop]` — Closed loops, solid on the left: outer loops counter-clockwise, holes clockwise (seen from the normal).
- `.normal: Vec3` — Plane normal.
- `.origin: Vec3` — Origin of the plane coordinates.
- `.s_axis: Vec3` — First in-plane unit axis s.
- `.t_axis: Vec3` — Second in-plane unit axis t = normal x s.

### `Plane`

Infinite plane through `point` with normal `normal`.

- `Plane(point: Vec3, normal: Vec3 = (0, 0, 1))` — Plane through `point` with the given `normal` (normalised).
- `.distance(p: Vec3) -> float` — Signed distance of p from the plane (positive on the side of the normal).
- `.normal: Vec3` (settable) — Unit normal.
- `.point: Vec3` (settable) — A point on the plane (origin of its in-plane coordinates).

### `PlaneFit`

Best-fit plane (orthogonal regression).

- `.max_deviation: float` — largest distance of a point from the plane
- `.plane: Plane` — plane through the centroid, unit normal
- `.rms: float` — root mean square of the distances

### `ProjectedCurve`

A curve on a surface: 3D curve and p-curve.

- `.curve: NurbsCurve` — 3D curve on the surface.
- `.pcurve: NurbsCurve` — (u, v, 0) in the surface's parameter domain (same parameter as curve).
- `__iter__() -> Iterator` — Iteration `for x in a`.

### `RebuiltCurve`

Result of `rebuild_curve`: the rebuilt curve and its max error.

- `.curve: NurbsCurve` — the rebuilt non-rational B-spline curve
- `.max_error: float` — max distance to the input at equal parameters (dense samples)

### `RebuiltSurface`

Result of `rebuild_surface`: the rebuilt surface and its max error.

- `.max_error: float` — max distance to the input at equal parameters (dense grid)
- `.surface: NurbsSurface` — the rebuilt non-rational B-spline surface

### `SectionCurve`

Piece of a section loop on one face of a shell.

- `.curve: NurbsCurve` — In plane coordinates (s, t, 0).
- `.curve3d: NurbsCurve` — The same piece in 3D (same parameter as `curve`).
- `.exact: bool` — True for exact isoparametric curves of the face; otherwise traced and fitted.
- `.face: int` — Index of the source face in the shell.
- `.normal(t: float) -> Vec3` — Outward unit normal of the shell at parameter t.
- `.pcurve: NurbsCurve` — (u, v, 0) on the source face's surface.

### `SectionLoop`

One closed loop of a planar section: consecutive curves and the signed enclosed area.

- `.area: float` — Signed area in (s, t): > 0 outer loop, < 0 hole.
- `.curves: list[SectionCurve]` — Consecutive `SectionCurve` pieces.

### `SectionOptions`

Options of `section`.

- `SectionOptions()` — Default options; set the attributes afterwards.
- `.fit_tolerance: float` (settable) — Accuracy of traced and fitted section curves on general faces (model units); 0 = 1e-9 of the model size.
- `.join_tolerance: float` (settable) — Chaining tolerance (model units, default 1e-6): ends of neighbouring pieces closer than this are joined.
- `.s_axis: Vec3` (settable) — First in-plane axis s (projected into the plane and normalised); zero vector = automatic.

### `Shell`

Trimmed faces sewn along shared edges.

- `.area(rel_tol: float = 1e-08) -> float` — Total area of all faces (adaptive quadrature to `rel_tol`).
- `.components() -> list[Shell]` — Connected components as separate shells (outer boundary, void boundaries, lumps).
- `.edges: list[ShellEdge]` — Edges: each `ShellEdge` has its 3D curve, its two vertex indices and its uses (face, loop, curve, direction).
- `.face_labels: list[str]` — The label of every face, in face order.
- `.faces: list[TrimmedSurface]` — The trimmed faces.
- `.faces_labeled(label: str) -> list[int]` — Indices of the faces with this label.
- `.free_edges() -> list[int]` — Indices of edges used by only one face (boundary of an open shell).
- `.is_closed: bool` — True if every edge is used exactly twice (watertight).
- `.is_consistently_oriented: bool` — True if every edge is used in opposite directions by its two faces.
- `.mass_properties(rel_tol: float = 1e-08) -> MassProperties` — Volume, centroid and inertia tensor (about the centroid, unit density) of the solid.
- `.non_manifold_edges() -> list[int]` — Indices of edges used by more than two faces.
- `.oriented() -> Shell` — The shell with faces flipped so that the orientation is consistent and every closed component has outward normals (voids inward).
- static `Shell.sew(faces: Sequence[TrimmedSurface], tolerance: float = 1e-06) -> Shell` — Sew trimmed faces into a shell: boundary curves closer than `tolerance` become shared edges (split at T-junctions); vertices are merged.
- `.transformed(transform: Transform) -> Shell` — The shell mapped by an affine transform (mirrors keep the normals outward).
- `.vertices: list[Vec3]` — Vertex positions.
- `.volume(rel_tol: float = 1e-08) -> float` — Enclosed volume (divergence theorem over the faces, adaptive to `rel_tol`); requires a closed shell. Positive for outward normals.
- `.with_face_labels(labels: Sequence[str]) -> Shell` — The same shell with one label per face (they survive Booleans, transforms and STEP).
- `__rmul__(arg0: Transform) -> Shell` — Operator `b * a`.

### `SilhouetteCurve`

Silhouette curve of a surface for a view; curve and pcurve share their parameter.

- `.closed: bool` — True if the 3D curve is closed.
- `.crease: bool` — crease of the surface between a part facing the viewer and one facing away
- `.curve: NurbsCurve` — 3D curve on the surface
- `.cusps: list[float]` — parameters where the view ray is tangent to the curve (cusps of the projection)
- `.face: int` — index of the source face (shells)
- `.pcurve: NurbsCurve` — (u, v, 0) on the surface

### `SurfaceChamfer`

Chamfer between two surfaces.

- `.contact1: NurbsCurve` — Contact curve on s1 (the boundary u = 0 of the chamfer).
- `.contact2: NurbsCurve` — Contact curve on s2 (the boundary u = 1 of the chamfer).
- `.edge: NurbsCurve` — the intersection curve that is cut off
- `.pcurve1: NurbsCurve` — contact1 as (u, v, 0) on s1
- `.pcurve2: NurbsCurve` — contact2 as (u, v, 0) on s2
- `.surface: NurbsSurface` — u in [0, 1]: straight lines from the contact on s1 to the contact on s2; v along the edge

### `SurfaceContinuity`

Continuity of two surfaces along a common edge.

- `.geometric: int` — -1 no contact, 0 position (G0), 1 tangent plane (G1), 2 curvature (G2)
- `.max_angle: float` — radians
- `.max_curvature_difference: float` — Largest difference (Frobenius norm) of the curvature tensors of the two surfaces along the edge.
- `.max_gap: float` — Largest distance between the two surfaces at the samples of the edge.

### `SurfaceCurvature`

Curvature of a surface at a point: Gaussian and mean curvature, principal curvatures and directions, normal.

- `.dir1: Vec3` — Unit principal direction of k1.
- `.dir2: Vec3` — Unit principal direction of k2.
- `.gaussian: float` — Gaussian curvature K = k1 k2.
- `.k1: float` — Larger principal curvature.
- `.k2: float` — Smaller principal curvature.
- `.mean: float` — Mean curvature H = (k1 + k2) / 2 (positive where the surface bends towards the normal).
- `.normal: Vec3` — Unit normal.

### `SurfaceDirection`

Parameter directions and orientation of a surface.

- `.closed_u: bool` — True if the sides u_min and u_max coincide.
- `.closed_v: bool` — True if the sides v_min and v_max coincide.
- `.normal: Vec3` — Unit normal S_u × S_v at the middle of the domain.
- `.planar: bool` — True if the control points lie in a plane within the tolerance.
- `.point: Vec3` — Surface point at the middle of the parameter domain.
- `.poles: tuple[bool, bool, bool, bool]` — sides collapsing to a point: (u_min, u_max, v_min, v_max)
- `.u_dir: Vec3` — Unit direction of increasing u there.
- `.v_dir: Vec3` — Unit direction of increasing v there.

### `SurfaceFillet`

Rolling-ball fillet between two surfaces.

- `.contact1: NurbsCurve` — Contact curve on the first surface (= fillet boundary u = 0).
- `.contact2: NurbsCurve` — Contact curve on the second surface (= fillet boundary u = 1).
- `.pcurve1: NurbsCurve` — contact1 as (u, v, 0) on s1
- `.pcurve2: NurbsCurve` — contact2 as (u, v, 0) on s2
- `.spine: NurbsCurve` — path of the ball centre
- `.surface: NurbsSurface` — u: arc from contact1 (u=0) to contact2 (u=1), v: spine parameter

### `SurfaceIntersection`

Result of a surface/surface intersection or a planar section: intersection curves and isolated touching points.

- `.curves: list[IntersectionCurve]` — Intersection branches (`IntersectionCurve`); closed surfaces split them at their seams.
- `.points: list[SurfaceSurfacePoint]` — Isolated touching points.

### `SurfaceIntersectionOptions`

Options of the surface/surface intersection.

- `SurfaceIntersectionOptions(tolerance: Tolerance = Tolerance(), chord_tolerance: float = 0.0, fit_tolerance: float = 0.0)` — 0 for `chord_tolerance` / `fit_tolerance` selects automatic values from the model size.
- `.chord_tolerance: float` (settable) — Maximal deviation of the traced polyline from the true curve (step control); 0 = automatic.
- `.fit_tolerance: float` (settable) — Accuracy of the fitted 3D curve and p-curves; 0 = automatic.
- `.tolerance: Tolerance` (settable) — Basic tolerances for point convergence.

### `SurfaceMapRange`

Smallest and largest sample of a map with their parameters.

- `.max: float` — Largest sample.
- `.max_uv: tuple[float, float]` — Parameters (u, v) of the largest sample.
- `.min: float` — Smallest sample.
- `.min_uv: tuple[float, float]` — Parameters (u, v) of the smallest sample.

### `SurfacePoint`

Result of a closest-point query on a surface.

- `.distance: float` — Distance from the query point.
- `.point: Vec3` — The closest point.
- `.u: float` — u parameter of the closest point.
- `.v: float` — v parameter of the closest point.

### `SurfaceSurfacePoint`

A point on both surfaces with its parameters on each.

- `.point: Vec3` — The 3D point.
- `.u1: float` — u on the first surface.
- `.u2: float` — u on the second surface.
- `.v1: float` — v on the first surface.
- `.v2: float` — v on the second surface.

### `SweepOptions`

Options of `sweep`, `sweep2` and `rail_revolve`.

- `SweepOptions(frame: SweepFrame = SweepFrame.RotationMinimizing, twist: float = 0.0, tolerance: float = 1e-06, max_sections: int = 1025)` — All arguments are keywords.
- `.frame: SweepFrame` (settable) — Orientation of the profile along the trajectory (`SweepFrame`).
- `.max_sections: int` (settable) — Upper limit for the number of sections inserted while refining.
- `.tolerance: float` (settable) — Deviation of the swept surface from the ideal sweep, verified between the sections.
- `.twist: float` (settable) — Additional rotation about the tangent at the end of the path (radians), distributed linearly along the arc length.

### `TangentLine`

Line segment tangent to curves.

- `.line: NurbsCurve` — segment from point_a to point_b
- `.param_a: float` — tangent parameter on the first curve (nan if it is a point)
- `.param_b: float` — tangent parameter on the second curve
- `.point_a: Vec3` — line start: tangent point on the first curve, or the given point
- `.point_b: Vec3` — line end: tangent point on the second curve

### `Tolerance`

Per-operation tolerances (length, parametric, angular). Tolerances of iterative algorithms: model-space length, parameter differences and angles. The defaults are the kernel constants (1e-10, 1e-10, 1e-9).

- `Tolerance(length: float = 1e-10, parametric: float = 1e-10, angular: float = 1e-09)` — All arguments are keywords.
- `.angular: float` (settable) — Angular tolerance (radians).
- `.length: float` (settable) — Model-space length tolerance (model units).
- `.parametric: float` (settable) — Parameter-space tolerance.

### `Transform`

Affine transformation (3 x 4 matrix, row-major). Compose with `t1 * t2` (t2 applied first); apply to geometry with `obj.transformed(t)` or `t * shell`.

- `Transform()`  
  `Transform(matrix_3x4_row_major: Sequence[float] (12))` — Identity, or from the 12 entries of a row-major 3 x 4 matrix [R | t].
- `.apply(arg0: Vec3) -> Vec3` — Apply to a point (rotation/scale and translation).
- `.apply_vector(arg0: Vec3) -> Vec3` — Apply to a direction (linear part only, no translation).
- `.matrix() -> list[float] (12)` — The 12 entries of the row-major 3 x 4 matrix.
- static `Transform.mirror(point: Vec3, normal: Vec3) -> Transform` — Reflection in the plane through `point` with the given `normal`.
- static `Transform.rotation(origin: Vec3, axis: Vec3, angle: float) -> Transform` — Rotation by `angle` radians about the axis through `origin` along `axis` (right-hand rule).
- static `Transform.scaling(center: Vec3, factor: float) -> Transform`  
  `Transform.scaling(center: Vec3, factors: Vec3) -> Transform` — Uniform (`factor`) or per-axis (`factors`) scaling about `center`.
- static `Transform.translation(d: Vec3) -> Transform` — Translation by d.
- `__mul__(arg0: Transform) -> Transform` — Operator `a * b`.

### `TrimLoop`

Closed loop of parameter-space curves (x = u, y = v, z = 0). Closed loop of curves in the (u, v) parameter space of a surface (curves with z = 0, consecutive end points coincide). Outer loops run counter-clockwise, holes clockwise.

- `TrimLoop(curves: Sequence[NurbsCurve])` — Create a loop from consecutive curves (u, v, 0); it must be closed.
- `.curves: list[NurbsCurve]` — The curves of the loop.
- `.polyline(samples_per_span: int = 16) -> list[list[float] (2)]` — Sampled (u, v) points of the loop, `samples_per_span` per knot span.
- `.reversed() -> TrimLoop` — The loop traversed in the opposite direction.
- `.signed_area() -> float` — Signed area enclosed in parameter space, from the polyline (approximate; positive for counter-clockwise loops). For exact areas use `TrimmedSurface.area` or `region_area`.

### `TrimmedSurface`

A NURBS surface restricted to the region inside an outer trim loop and outside its holes (all in parameter space). Faces of shells are trimmed surfaces.

- `TrimmedSurface(surface: NurbsSurface)`  
  `TrimmedSurface(surface: NurbsSurface, outer: TrimLoop, holes: Sequence[TrimLoop] = [], tolerance: Tolerance = Tolerance())` — Untrimmed surface (the whole domain), or the surface with an outer loop and holes; the loops are checked against `tolerance`.
- `.area(rel_tol: float = 1e-09) -> float` — Area of the trimmed region (Green's theorem on the exact trim curves, adaptive to `rel_tol`).
- `.boundary_curves(tolerance: float = 1e-07) -> list[list[NurbsCurve]]` — The trim loops mapped to 3D: one list of curves per loop (outer first), exact where possible, otherwise fitted within `tolerance`.
- `.contains(u: float, v: float) -> bool` — True if the parameter point (u, v) lies in the trimmed region.
- `.flipped() -> TrimmedSurface` — The same face with opposite orientation (normal flipped).
- `.holes: list[TrimLoop]` — The hole loops.
- `.is_untrimmed: bool` — True if the face covers the whole parameter domain.
- `.label: str` — Face label (boundary condition / physical group); kept by Booleans, transforms, STEP.
- `.outer: TrimLoop` — The outer loop.
- `.surface: NurbsSurface` — The underlying surface.
- `.transformed(transform: Transform) -> TrimmedSurface` — The face mapped by an affine transform (loops and label kept; reflections keep the normal by flipping the face).
- `.with_hole(hole: TrimLoop) -> TrimmedSurface` — A copy with one more hole loop.
- `.with_label(label: str) -> TrimmedSurface` — The same face with the given label.
- `__rmul__(arg0: Transform) -> TrimmedSurface` — Operator `b * a`.

### `UnrolledSurface`

Result of unroll().

- `.max_error: float` — Largest deviation of the pattern from the exact development at the check points of the fit.
- `.pattern: NurbsSurface` — flat pattern in the plane z = 0; pattern(u, v) is the image of S(u, v)

### `Vec3`

3D point or vector with components `x`, `y`, `z`; supports `+`, `-`, scalar `*`, indexing and iteration. Tuples and lists of three numbers convert implicitly wherever a `Vec3` is expected.

- `Vec3()`  
  `Vec3(x: float, y: float, z: float)`  
  `Vec3(arg0: Sequence[float] (3))` — Zero vector, from three components, or from a sequence of three numbers.
- `.to_tuple() -> tuple[float, float, float]` — The components as a tuple (x, y, z).
- `.x: float` (settable) — x component.
- `.y: float` (settable) — y component.
- `.z: float` (settable) — z component.
- `__add__(arg0: Vec3) -> Vec3` — Operator `a + b`.
- `__getitem__(arg0: int) -> float` — Indexing `a[i]`.
- `__iter__() -> Iterator` — Iteration `for x in a`.
- `__len__() -> int` — `len(a)`.
- `__mul__(arg0: float) -> Vec3` — Operator `a * b`.
- `__neg__() -> Vec3` — Operator `-a`.
- `__rmul__(arg0: float) -> Vec3` — Operator `b * a`.
- `__sub__(arg0: Vec3) -> Vec3` — Operator `a - b`.

### `XYPolynomialSag`

Freeform sag: base asphere plus \sum c_{ij} (x/r_n)^i (y/r_n)^j with the normalisation radius r_n.

- `XYPolynomialSag(base: Asphere = Asphere(), norm_radius: float = 1.0, terms: Sequence[tuple[int, int, float]] = [])` — terms: list of (i, j, c) for c * (x/norm)^i * (y/norm)^j
- `.sag(x: float, y: float) -> float` — Sag z at (x, y).

### `ZernikeSag`

Freeform sag: base asphere plus \sum_j c_j Z_j(\rho, \theta) with \rho = r / r_n.

- `ZernikeSag(base: Asphere = Asphere(), norm_radius: float = 1.0, indexing: ZernikeIndexing = ZernikeIndexing.Fringe, coeffs: Sequence[float] = [])` — `coeffs[0]` is the coefficient of Z_1 in the chosen indexing.
- `.base: Asphere` (settable) — Base asphere.
- `.coeffs: list[float]` (settable) — Zernike coefficients c_1, c_2, …
- `.indexing: ZernikeIndexing` (settable) — Ordering and normalisation of the terms (`ZernikeIndexing`).
- `.norm_radius: float` (settable) — Normalisation radius r_n.
- `.sag(x: float, y: float) -> float` — Sag z at (x, y).

### Enumerations

- **`BooleanOp`** — Boolean operation for `boolean`. Values: `Union` (Material of a or b.); `Intersection` (Material of a and b.); `Difference` (Material of a without b.)
- **`CurveEnd`** — End of a curve for `extend`. Values: `Start` (At the start of the domain.); `End` (At the end of the domain.)
- **`DrawingCurveKind`** — Source of a line of a drawing. Values: `Edge` (Edge of a shell.); `Silhouette` (Silhouette curve of a face.)
- **`Extension`** — How `extend` continues the geometry. Values: `Natural` (Polynomial (rational) continuation of the boundary B\'ezier piece: C^∞ at the old boundary, conics stay conics; fails if weights would become ≤ 0.); `Linear` (Straight continuation along the boundary tangent: G1, and C1 for non-rational geometry.)
- **`Parametrization`** — Choice of the data parameters of fitting and lofting (P&T Eq. 9.3–9.6). Values: `Uniform` (Equally spaced.); `ChordLength` (Proportional to the distance between consecutive points.); `Centripetal` (Proportional to the square root of the distance (robust at sharp turns; the default).)
- **`RulingDirection`** — Parameter direction of the straight lines of a ruled surface. Values: `None_` (Not ruled along a parameter direction.); `U` (The lines v = const are straight (the surface is linear in u).); `V` (The lines u = const are straight (the surface is linear in v).)
- **`SurfaceSide`** — Boundary of a surface for `extend`. Values: `UMin` (Boundary u = u_min.); `UMax` (Boundary u = u_max.); `VMin` (Boundary v = v_min.); `VMax` (Boundary v = v_max.)
- **`Sweep2Scaling`** — Scaling of the sections of `sweep2`. Values: `Uniform` (The whole section is scaled by the ratio of the rail distances.); `NonUniform` (Only the rail-to-rail direction is scaled; heights above it are kept.)
- **`SweepFrame`** — Orientation of the profile along the sweep path. Values: `RotationMinimizing` (Rotation-minimising frame (double reflection): no spurious twist, defined at inflections and on straight parts.); `Fixed` (Pure translation of the profile (no rotation).)
- **`ZernikeIndexing`** — Ordering and normalisation of Zernike polynomials. Values: `Fringe` (Fringe / University of Arizona ordering, unnormalised (Z_4 = 2\rho^2 - 1).); `Noll` (Noll ordering, normalised to unit RMS over the unit disk.)

## `nurbforge.mesh`

Meshing via gmsh (optional dependency: `pip install nurbforge[mesh]`).

### `mesh.GmshSession`

Context manager around `gmsh.initialize()` / `gmsh.finalize()`. Shells added with `add_shell` keep their face labels, which `face_tags` maps to gmsh surfaces and `refine` uses for local mesh sizes.

- `GmshSession(model_name: str = nurbforge, verbose: bool = False)` — Initialize self. See help(type(self)) for accurate signature.
- `.add_curve(curve: NurbsCurve) -> int` — Add a curve exactly (`add_curve`); returns the gmsh curve tag.
- `.add_shell(shell: Shell, solid: bool = True, physical_groups: bool = False, volume_name: str = volume) -> list[int]` — As the module function `add_shell`; in addition the session remembers the face labels for `face_tags` and `refine`.
- `.add_surface(surface: NurbsSurface | TrimmedSurface) -> int` — Add a surface or trimmed surface exactly; returns the gmsh surface tag.
- `.face_tags() -> dict[str, list[int]]` — The current gmsh surface tags of every face label of the shells added with `add_shell` (unlabelled faces: `"unlabeled"`).
- `.mesh(dim: int = 2, size: float | None = None, order: int = 1, heal: bool = True, curvature: int | None = None) -> None` — Synchronise the OCC model and generate a mesh. `size` is the global maximum element size, `curvature` the number of elements per 2*pi of curvature (smaller elements on strongly curved faces), local sizes come from `refine`.
- `.refine(where: str | int | Iterable[str | int], size: float, distance: float = 0.0) -> None` — Local mesh size `size` on the given faces: face labels of shells added with `add_shell` and/or gmsh surface tags. With `distance > 0` the size grows linearly from `size` on the faces to the global size (`mesh(size=...)`) at that distance, also into the volume; with `distance = 0` only the faces themselves are refined. Applied by `mesh`; several refinements combine by taking the smallest size.
- `.surface_mesh(orient: bool = True) -> SurfaceMesh` — The triangle surface mesh of the model as a `SurfaceMesh`, for boundary element solvers. Raises `MeshTopologyError` unless the mesh is watertight (every edge shared by exactly two triangles) and orientable. With `orient` the triangles are made consistent across surfaces and their normals point out of the solid (into the void for a cavity); the gmsh mesh is reversed accordingly, so `write` exports the same orientation. Call after `mesh(dim=2)`; the triangle groups follow the face labels.
- `.write(path: str) -> None` — Write mesh/geometry; format from extension (.msh, .vtk, .stl, .step, .brep ...).

- `nurbforge.mesh.bem_mesh(shell: Shell, size: float | None = None, order: int = 1, curvature: int | None = None, refine: dict | None = None, path: str | None = None) -> SurfaceMesh` — One-shot helper for boundary element solvers: mesh the closed `shell` with triangles and return a watertight, outward-oriented `SurfaceMesh` (see `GmshSession.surface_mesh`). `refine` maps face labels or surface tags to local sizes; `path` additionally writes the mesh (.msh, .stl, .vtk, ...) with physical groups per face label.
### `mesh.SurfaceMesh`

A triangle surface mesh with 0-based connectivity. `nodes` (N, 3) coordinates; `triangles` (M, 3) for linear or (M, 6) for quadratic elements (gmsh node order: corners, then the midpoints of edges 0-1, 1-2, 2-0); `surfaces` (M,) the gmsh surface tag of every triangle; `node_tags` (N,) the gmsh node tags; `groups` triangle indices per face label (shells added with `add_shell`).

- `SurfaceMesh(nodes: np.ndarray, triangles: np.ndarray, surfaces: np.ndarray, node_tags: np.ndarray, groups: dict[str, np.ndarray] = <factory>)` — Initialize self. See help(type(self)) for accurate signature.
- `.areas() -> np.ndarray` — Areas of the (flat) corner triangles, (M,).
- `.boundary_edges() -> np.ndarray` — Corner-node edges used by only one triangle, (K, 2); empty for a watertight mesh.
- `.centroids() -> np.ndarray` — Centroids of the corner triangles, (M, 3).
- `.is_consistently_oriented() -> bool` — Every directed edge occurs once, i.e. neighbours traverse a shared edge in opposite directions.
- `.is_watertight() -> bool` — Every edge is shared by exactly two triangles.
- `.normals() -> np.ndarray` — Unit normals of the (flat) corner triangles, (M, 3).
- `.order` — Element order: 1 for 3-node, 2 for 6-node triangles.
- `.volume() -> float` — Enclosed volume (divergence theorem over the corner triangles); positive for an outward-oriented closed mesh, cavities count negative.
- `__eq__(other)` — Comparison `a == b`.

### `mesh.MeshTopologyError` (exception, subclass of ValueError)

The surface mesh is not a closed, orientable 2-manifold.

- `nurbforge.mesh.add_curve(curve: NurbsCurve) -> int` — Add a NURBS curve to the current gmsh OCC model. Returns the curve tag.
- `nurbforge.mesh.add_surface(surface: NurbsSurface, wire_tags: list[int] | None = None) -> int` — Add a NURBS surface to the current gmsh OCC model. Returns the surface tag. gmsh expects control points with the u index running fastest, while Nurbforge stores them row-major with v fastest, hence the reordering below. `wire_tags` (outer contour first, then holes) trims the surface with wires defined in its parameter space.
- `nurbforge.mesh.add_trimmed_surface(trimmed: TrimmedSurface) -> int` — Add a trimmed NURBS surface: the trim loops are passed exactly as B-spline wires in the parameter space of the surface (outer contour first, then holes). Returns the tag.
- `nurbforge.mesh.add_shell(shell: Shell, solid: bool = True, physical_groups: bool = False, volume_name: str = volume) -> list[int]` — Add a shell with its B-rep topology: faces share their edges. For a closed shell and `solid=True` a volume is created. Returns the gmsh surface tag of every face (in the order of `shell.faces`, without repetitions) and the volume tag last, if one was created. The shell is handed over through Nurbforge's native STEP writer (exact NURBS, shared edges) and OCC's STEP reader, so OCC builds the topology itself – including the seam edges of closed surfaces (cylinders, cones, spheres, tori), which a face-by-face transfer with duplicate removal or sewing cannot reproduce once such a face is trimmed (e.g. a bore or a slot from a Boolean). `physical_groups=True` creates one physical surface group per face label (`TrimmedSurface.label`, e.g. set with `Shell.with_face_labels` before a Boolean; unlabelled faces go to `"unlabeled"`) and a physical volume `volume_name`, so boundary conditions can be assigned by name in the solver. Faces are matched to the gmsh surfaces by an interior point.
- `nurbforge.mesh.export_step(objects, path: str, solid: bool = True) -> None` — Write surfaces, trimmed surfaces and shells to a STEP (or BREP/IGES, by extension) file via gmsh's OpenCASCADE kernel. Closed shells become solids if `solid` is set.
- `nurbforge.mesh.mesh_surfaces(surfaces: Iterable[NurbsSurface | TrimmedSurface], path: str, size: float | None = None, order: int = 1) -> None` — One-shot helper: add surfaces, mesh in 2D and write to `path`.

## `nurbforge.llm`

LLM interface (Ollama) for natural-language geometry construction. Design principle: the model never executes arbitrary code. It can only call the whitelisted tools in `nurbforge.llm.tools`, which operate on a named object store (`nurbforge.llm.tools.Workspace`).

### `llm.GeometryAgent`

Tool-calling agent. After the model's final answer, up to `verify_rounds` verification rounds follow: the objects created or changed for the request are measured (`verify.report`: bounding boxes, lengths, areas, volumes, closedness, detected issues such as open shells) and the report is sent back to the model, which either fixes the geometry with further tool calls or confirms. The rounds stop as soon as the model confirms without calling tools. The last report is kept in `verification`.

- `GeometryAgent(client: OllamaClient | None = None, workspace: Workspace | None = None, max_steps: int = 20, verify_rounds: int = 1)` — Initialize self. See help(type(self)) for accurate signature.
- `.ask(prompt: str) -> str` — Send a request to the model, run the tool calls it makes (and the verification rounds) and return its final answer.

### `llm.OllamaClient`

Thin wrapper around the Ollama REST API (`/api/chat`, `/api/tags`).

- `OllamaClient(host: str = http://localhost:11434, model: str = qwen2.5:14b, timeout: float = 300.0)` — Initialize self. See help(type(self)) for accurate signature.
- `.chat(messages: list[dict[str, Any]], tools: list[dict[str, Any]] | None = None, temperature: float = 0.0) -> dict[str, Any]` — Single non-streaming chat turn. Returns the assistant `message` dict.
- `.list_models() -> list[str]` — Names of the models available on the Ollama server.

### `llm.OllamaError` (exception, subclass of RuntimeError)

Raised when the Ollama server cannot be reached or returns an error.

### `llm.Workspace`

Named object store the LLM operates on.

- `Workspace(objects: dict[str, Any] = <factory>)` — Initialize self. See help(type(self)) for accurate signature.
- `.get(name: str) -> Any` — The object stored under `name`; raises `KeyError` listing the known names.
- `.put(name: str, obj: Any) -> str` — Store `obj` under `name` (replacing an existing one) and return a short confirmation for the model.
- `__eq__(other)` — Comparison `a == b`.

- `nurbforge.llm.tools.call_tool(ws: Workspace, name: str, arguments: dict[str, Any]) -> str` — Execute a tool call; errors are returned as text so the LLM can correct itself.
- `nurbforge.llm.verify.inspect_object(obj: Any) -> dict[str, Any]` — Measured properties of a workspace object (bounding boxes from dense samples).
- `nurbforge.llm.verify.issues(props: dict[str, Any]) -> list[str]` — Sanity problems of an object, from its measured properties.
- `nurbforge.llm.verify.describe(name: str, obj: Any) -> str` — One line per object: measured properties and issues.
- `nurbforge.llm.verify.report(objects: dict[str, Any]) -> str` — Verification report for several named objects.
- `nurbforge.llm.verify.check(obj: Any, prop: str, expected: Any, tolerance: float | None = None) -> str` — Compare a measured property (e.g. 'volume', 'bbox.max', 'closed') with an expectation. Numbers (and number lists) pass within `tolerance` (default: 1e-6 relative, at least 1e-9 absolute); other values must be equal.
- `nurbforge.llm.constraints.evaluate(expr: str, variables: dict[str, float]) -> float` — Evaluates an arithmetic expression over the variables (whitelist, no Python code).
- `nurbforge.llm.constraints.substitute(value: Any, variables: dict[str, float]) -> Any` — Replaces "$name" and "${expression}" strings (also inside lists / dicts) by numbers.
- `nurbforge.llm.constraints.solve(ws: Any, recipe: list[dict[str, Any]], unknowns: dict[str, dict[str, float]], targets: list[dict[str, Any]], call_tool, tolerance: float = 1e-08, max_iterations: int = 40) -> dict[str, Any]` — Bounded Levenberg–Marquardt on the relative target residuals; on success the recipe is executed in `ws` with the solution. Returns {"converged", "values", "residuals", "iterations"}.

## LLM tools (`nurbforge.llm.tools.TOOLS`)

JSON tools for an LLM agent working on a `Workspace` of named objects; call them with `nurbforge.llm.tools.call_tool(ws, name, arguments)`. Angles in degrees. `name?` = optional parameter.

- `make_line(name: string, p0: number[3], p1: number[3])` — Create a straight line segment.
- `make_circle(name: string, center: number[3], normal: number[3], radius: number)` — Create a full circle (exact NURBS).
- `make_arc(name: string, center: number[3], x_axis: number[3], y_axis: number[3], radius: number, start_deg: number, end_deg: number)` — Create a circular arc. Angles in degrees, measured from x_axis to y_axis.
- `make_polygon(name: string, center: number[3], normal: number[3], radius: number, n: integer)` — Create a closed regular polygon (n vertices on a circle of the given radius about center, in the plane with the given normal).
- `make_helix(name: string, axis_point: number[3], axis_dir: number[3], start: number[3], pitch: number, turns: number, end_radius?: number, left_handed?: boolean)` — Create a helix about an axis starting at a point (its distance from the axis is the radius): pitch = rise per turn (0: planar spiral), turns > 0, optional end_radius for a conical helix or spiral, left_handed to wind clockwise.
- `round_corners(name: string, source: string, radius: number)` — Round all corners of a planar curve (polyline, polygon, rectangle outline) with circular arcs of the given radius; one smooth curve.
- `chamfer_corners(name: string, source: string, distance: number)` — Chamfer all corners of a curve (polyline, polygon, rectangle outline) with straight cuts at the given distance from each corner; one curve.
- `make_rectangle(name: string, origin: number[3], edge_u: number[3], edge_v: number[3])` — Create a planar rectangular surface.
- `revolve(name: string, profile: string, axis_point: number[3], axis_dir: number[3], angle_deg: number)` — Revolve a curve about an axis to create a surface. Angle in degrees.
- `extrude(name: string, profile: string, direction: number[3])` — Extrude a curve along a direction vector to create a surface.
- `ruled(name: string, curve0: string, curve1: string)` — Create a ruled surface between two curves.
- `translate(name: string, source: string, vector: number[3])` — Copy an object translated by a vector.
- `rotate(name: string, source: string, axis_point: number[3], axis_dir: number[3], angle_deg: number)` — Copy an object rotated about an axis. Angle in degrees.
- `mirror(name: string, source: string, plane_point: number[3], plane_normal: number[3])` — Copy an object mirrored at a plane (point + normal).
- `polar_array(name: string, source: string, axis_point: number[3], axis_dir: number[3], count: integer, total_angle_deg: number)` — Create count rotated copies about an axis; stored as name_0..name_{count-1}.
- `linear_array(name: string, source: string, step: number[3], count: integer)` — Create count translated copies; stored as name_0..name_{count-1}.
- `list_objects()` — List all objects in the workspace.
- `closest_point(name: string, point: number[3])` — Closest point of an object to a given point, with its parameters and the distance.
- `interpolate_points(name: string, points: number[3][], degree?: integer)` — Create a smooth curve through the given points (in order).
- `fit_surface(name: string, points: number[3][], smoothing?: number, num_ctrl?: integer)` — Fit a smooth surface to unordered measured points (e.g. a height map). smoothing: 0.001 follows the points closely, 0.1-1 smooths noisy data. Reports the max and rms deviation.
- `loft(name: string, sections: string[], degree?: integer)` — Create a smooth surface through several section curves (in the given order).
- `sweep(name: string, profile?: string, path: string, twist_deg?: number, profiles?: string[])` — Sweep a profile curve along a path curve. The profile must be drawn at the start of the path (typically in the plane normal to the path). Optional twist in degrees over the whole path. With 'profiles' (several curves drawn where they cross the path, e.g. a circle at the start and a larger one at the end) the cross-section blends between them.
- `rail_revolve(name: string, profile: string, rail: string, axis_point: number[3], axis_dir: number[3], scale_height?: boolean)` — Revolve a profile curve about an axis while its outer end follows a rail curve (e.g. a reflector with an elliptic or square-ish rim): the profile is scaled radially so its end stays on the rail; scale_height also scales heights.
- `extrude_tapered(name: string, profile: string, direction: number[3], angle_deg: number)` — Tapered (draft-angle) extrusion of a planar profile along a direction normal to its plane (the vector length is the height). angle_deg > 0 widens the section towards the top, < 0 narrows it. Corners are mitred; a profile with corners gives an open shell of faces (use 'cap' to close it).
- `ribbon(name: string, curve: string, width: number, plane_normal?: number[3])` — Planar strip of the given width along a curve (to the left seen from the plane normal's tip; negative width: right). Corners are mitred.
- `extrude_to_point(name: string, curve: string, apex: number[3])` — Cone-like surface from a curve to an apex point (exact for circles and arcs).
- `triangle_patch(name: string, c0: string, c1: string, c2: string)` — Three-sided surface filling three boundary curves that form a closed loop (any orientation).
- `rebuild(name: string, source: string, num_u: integer, num_v?: integer, degree?: integer)` — Rebuild a curve or surface with a given number of control points (uniform cubic B-spline by default), e.g. to simplify or unify it; reports the maximum deviation.
- `coons(name: string, c0: string, c1: string, d0: string, d1: string)` — Create a surface filling four boundary curves that form a closed loop (any orientation).
- `asphere(name: string, radius: number, conic?: number, coefficients?: number[], clear_radius: number)` — Create an optical lens surface (sag along +z, vertex at the origin): radius of curvature R (0 = flat; positive: centre at +z), conic constant k and optional even aspheric coefficients A4, A6, ... over the clear radius.
- `extend(name: string, source: string, where: start | end | umin | umax | vmin | vmax, length: number, mode?: natural | linear)` — Extend a curve at its start or end (or a surface across a side umin/umax/vmin/vmax) by a length. mode 'natural' continues the shape smoothly, 'linear' continues straight along the tangent.
- `offset(name: string, source: string, distance?: number, plane_normal?: number[3], corner_distances?: number[4])` — Offset a curve or surface by a signed distance. Surfaces move along their normal. Curves must be planar and move within their plane, to the left of the direction of travel seen from plane_normal (default [0, 0, 1]; e.g. inwards for a counter-clockwise circle). Fails if the distance reaches a radius of curvature. Surfaces only: corner_distances gives a varying distance (e.g. a layer with a thickness gradient) at the corners (u0,v0), (u1,v0), (u0,v1), (u1,v1), interpolated bilinearly; it replaces distance.
- `fillet(name: string, a: string, b: string, radius: number, side_a?: -1 | 0 | 1, side_b?: -1 | 0 | 1, plane_normal?: number[3])` — Round the corner between two objects with a constant radius. Two planar intersecting curves: circular arc tangent to both, stored as name, with the trimmed curves stored as name_a and name_b (a, arc, b form a smooth chain); the corner is chosen automatically unless side_a/side_b are given (+1 left / -1 right of the curve direction seen from plane_normal). Two surfaces: rolling-ball fillet surface stored as name, contact curves as name_contact1/name_contact2; side_a/side_b = +1 (default) put the ball on the normal side of a/b, -1 on the other side.
- `chamfer(name: string, a: string, b: string, distance_a: number, distance_b?: number, side_a?: -1 | 1, side_b?: -1 | 1)` — Straight cut between two intersecting curves: from the point of a at distance_a from their intersection to the point of b at distance_b (default distance_a). The cut is stored as name, the trimmed curves as name_a and name_b (a, cut, b form a chain). Two intersecting surfaces: the chamfer surface between the curves at these distances from the intersection curve is stored as name, the contact curves as name_contact1/name_contact2; side_a/side_b = +1 (default) put it on the normal side of a/b, -1 on the other side.
- `intersect(a: string, b: string, name?: string)` — Intersection of two objects: points for curve/curve and curve/surface; for surface/surface the intersection curves are stored as name_1, name_2, ...
- `trim_circle(name: string, source: string, center_u: number, center_v: number, radius: number)` — Restrict a surface to a circular region of its parameter domain (e.g. a circular aperture of an optical surface: for sag/graph surfaces u = x, v = y).
- `split(name: string, surface: string, cutter: string)` — Split a surface along its intersection with another surface; the pieces are stored as name_1, name_2, ...
- `section(name: string, surface: string, point: number[3], normal: number[3])` — Cut a surface with a plane (point and normal); the section curves are stored as name_1, name_2, ...
- `sew(name: string, faces: string[])` — Sew surfaces (and trimmed surfaces) along shared edges into a shell; reports whether it is closed and, if so, the enclosed volume.
- `solid(name: string, kind: box | sphere | cylinder | cone | torus, point: number[3], axis?: number[3], size?: number[3], radius?: number, radius2?: number, height?: number)` — Create a closed solid (shell with outward normals): kind 'box' (corner = 'point', edge lengths = 'size'), 'sphere' (center, radius), 'cylinder' (base center, axis, radius, height), 'cone' (base, axis, radius = base radius, radius2 = top radius, 0 for an apex, height) or 'torus' (center, axis, radius = major, radius2 = tube).
- `mass_properties(name: string)` — Volume, centroid and inertia tensor (about the centroid, unit density; principal moments and axes) of a closed shell.
- `boolean(name: string, a: string, b: string, op: union | intersection | difference)` — Boolean operation of two closed shells (from 'sew'): union, intersection or difference (a minus b). Reports the volume of the result.
- `blend(name: string, a: string, end_a: start | end, b: string, end_b: start | end, continuity?: 1 | 2)` — Smooth transition curve from an end of curve a to an end of curve b (G1 or G2 continuity). Ends are 'start' or 'end'.
- `blend_surfaces(name: string, a: string, side_a: umin | umax | vmin | vmax, b: string, side_b: umin | umax | vmin | vmax, continuity?: 0 | 1 | 2, tension?: number)` — Smooth transition surface between a side of surface a and a side of surface b (sides umin/umax/vmin/vmax): continuity 0 (ruled), 1 (tangent planes match) or 2 (also curvature). tension scales how far the blend bulges.
- `match_surface(name: string, surface: string, side: umin | umax | vmin | vmax, target: string, target_side: umin | umax | vmin | vmax, continuity?: 0 | 1 | 2)` — Modify a surface near one of its sides so that it joins a side of a target surface smoothly (continuity 0: edges meet, 1: tangent, 2: curvature). The opposite side stays unchanged. Stored under name.
- `project(name: string, curve: string, surface: string, direction?: number[3])` — Project a curve onto a surface, along the surface normal or along a direction. Stores the 3D curve as name and its parameter-space curve as name + '_uv'.
- `project_to_plane(name: string, source: string, point: number[3], normal: number[3], direction?: number[3])` — Flatten a curve or surface onto a plane (through point, with normal) by parallel projection along direction (default: the plane normal). Exact.
- `offset_on_surface(name: string, curve_uv: string, surface: string, distance: number)` — Offset a curve lying on a surface within the surface (geodesic distance; positive = left of the curve seen from the surface normal). The curve is given by its parameter-space curve (e.g. the name + '_uv' stored by project). Stores the 3D curve as name and its parameter-space curve as name + '_uv'.
- `tween(name: string, a: string, b: string, count: integer)` — Intermediate curves between two curves a and b (running in the same direction), stored as name_1 .. name_count.
- `split_solid(name: string, solid: string, point: number[3], normal: number[3], cut_label?: string)` — Split a closed solid by a plane into the part below (opposite to the normal, stored as name + '_below') and above (name + '_above'). The new faces in the plane get the label cut_label (default 'cut'), e.g. for a symmetry condition.
- `cap(name: string, shell: string)` — Close the planar openings of an open shell with planar faces (label 'cap').
- `thicken(name: string, shell: string, thickness: number, symmetric?: boolean)` — Give an open smooth shell or a surface a wall thickness along its normal (negative: other side; symmetric: half on each side). Side walls get label 'wall'.
- `hollow(name: string, solid: string, thickness: number, open_labels?: string[])` — Hollow a smooth closed solid to a wall of the given thickness: an inner void, or open at the faces whose labels are listed in open_labels.
- `pipe(name: string, path: string, radius: number, end_radius?: number)` — Solid pipe along a path curve with a circular cross-section of the given radius (end_radius: linear taper). Faces labelled 'wall', 'inlet', 'outlet'.
- `merge_faces(name: string, solid: string)` — Defeaturing: merge adjacent faces of a solid that lie on the same surface and carry the same label (e.g. coplanar pieces left by a union).
- `round_edges(name: string, solid: string, label_a: string, label_b: string, size: number, kind?: fillet | chamfer, sizes?: object[])` — Fillet (round) or chamfer the edges of a solid between the faces labelled label_a and label_b (e.g. 'top' and 'side' for the rim of a cylinder, 'top' and 'front' for a box edge); '*' matches any label ('top', '*': all edges of the top face; '*', '*': every edge). Convex edges lose material, concave edges gain it. Where selected edges meet at a corner they are mitred, continue tangentially, or (all three edges of a corner) get a spherical / triangular corner patch. Optional variable size: radius / distance at points on the edges. New faces are labelled 'fillet' or 'chamfer'.
- `remove_holes(name: string, solid: string, max_size: number)` — Defeaturing: fill the through and blind holes of a closed solid that are smaller than max_size across.
- `export_step(objects: string[], path: string, unit?: mm | cm | m | inch)` — Write shells (from 'sew' or 'boolean'; surfaces are sewn into a shell each) to a STEP file with the native writer. Closed shells become solids.
- `import_step(name: string, path: string, unit?: mm | cm | m | inch)` — Read the solids and surface models of a STEP file with the native reader. The shells are stored as name, name_1, name_2, ... (one per solid / shell).
- `solve(recipe: object[], unknowns: object, targets: object[])` — Find unknown dimensions from measured targets. 'recipe' is a list of tool calls {tool, args} whose numeric arguments may be "$name" or "${expression}" (+ - * / **, pi, sqrt, sin/cos/tan in degrees); 'unknowns' maps each name to {start, min, max}; 'targets' are {name, property, value} on measured properties (as in 'check': volume, area, length, bbox.max, ...). The recipe is re-run until the targets are met, then executed for real.
- `inspect(name: string)` — Measured properties of an object (type, bounding box, length / area / volume, closed, free edges, ...) and detected issues. Use it to verify results.
- `check(name: string, property: string, expected: any, tolerance?: number)` — Verify a measured property of an object against an expected value, e.g. property 'volume', 'area', 'length', 'closed', 'bbox.min', 'bbox.max', 'faces'. Returns PASS or FAIL with the measured value.
- `measure(name: string)` — Length of a curve or area of a surface.
- `evaluate(name: string, u: number, v?: number)` — Evaluate a curve at u or a surface at (u, v).
