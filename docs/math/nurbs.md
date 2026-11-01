# NURBS fundamentals and conventions

## Definitions

A NURBS curve of degree *p* with control points **P**ᵢ, weights *w*ᵢ > 0 and knot vector
U = {u₀,…,u_{n+p+1}}:

    C(u) = Σᵢ N_{i,p}(u) wᵢ Pᵢ / Σᵢ N_{i,p}(u) wᵢ ,   u ∈ [u_p, u_{n+1}]

A tensor-product surface of degrees (p, q):

    S(u,v) = Σᵢ Σⱼ N_{i,p}(u) N_{j,q}(v) w_{ij} P_{ij} / Σᵢ Σⱼ N_{i,p}(u) N_{j,q}(v) w_{ij}

Internally all evaluation is done in homogeneous coordinates **P**ʷ = (wP, w) and projected.

## Conventions

| Item | Convention |
|---|---|
| Knot vectors | Non-decreasing, size `n_ctrl + degree + 1`, interior multiplicity ≤ degree. Clamped vectors are typical but not required. |
| Default domain | `[0, 1]` for all constructors in `primitives` / `construct` |
| Surface net | Row-major `index = i * num_v + j`; i along u |
| Revolution | u = circumferential (degree 2, rational), v = profile |
| Extrusion / ruled | u = profile, v = sweep direction (degree 1) |
| Angles | Radians in C++/Python API; **degrees** in LLM tools |
| Circle orientation | Counter-clockwise about the given normal; start point on `any_perpendicular(normal)` |
| Tolerances | `nf::tol::kLength = 1e-10` (model units), `kParametric = 1e-10`, `kAngular = 1e-9` rad – see [ADR 0005](../adr/0005-tolerance-model.md) |

## Implemented algorithms (Piegl & Tiller numbering)

| Algorithm | Function |
|---|---|
| A2.1 FindSpan | `find_span` |
| A2.2 BasisFuns | `basis_functions` |
| A2.3 DersBasisFuns | `basis_function_derivatives` |
| A4.1 / A4.3 rational point | `NurbsCurve::evaluate`, `NurbsSurface::evaluate` |
| A4.2 RatCurveDerivs | `NurbsCurve::derivatives` |
| A3.6 SurfaceDerivsAlg1 + A4.4 RatSurfaceDerivs | `NurbsSurface::derivatives` (any order), `first_derivatives` |
| A5.1 CurveKnotIns | `NurbsCurve::insert_knot` |
| A5.4 RefineKnotVectCurve | `NurbsCurve::refine_knots` |
| A5.5 RefineKnotVectSurface | `NurbsSurface::refine_knots_u/v` (A5.4 per column/row) |
| A5.8 RemoveCurveKnot | `NurbsCurve::try_remove_knot` / `remove_knot`, `NurbsSurface::remove_knot_u/v` (common count over all rows) |
| Degree reduction (Eq. 5.41/5.42) | `NurbsCurve::reduce_degree`, `NurbsSurface::reduce_degree_u/v` – see below |
| A5.9 DegreeElevateCurve | `NurbsCurve::elevate_degree` (input clamped first) |
| A5.10 DegreeElevateSurface | `NurbsSurface::elevate_degree_u/v` (A5.9 per column/row) |
| A5.6 DecomposeCurve (equivalent) | `NurbsCurve::decompose_bezier`, `NurbsSurface::decompose_bezier` – refine every domain knot to multiplicity p, then slice; also works for unclamped knot vectors |
| A7.1 MakeNurbsCircle | `make_arc`, `make_circle` |
| A8.1 MakeRevolvedSurf | `revolve` |

## Knot removal and degree reduction tolerances

Tolerances of `remove_knot` and `reduce_degree` are **Cartesian deviations** (model units).
Internally they are checked on homogeneous control points: for non-rational geometry the
control-point deviation bounds the curve deviation directly; for rational geometry the tolerance
is scaled by w_min / (1 + max|P_i|) (P&T Eq. 5.30), which is conservative.

Degree reduction is not a port of A5.11. Instead:

1. decompose into Bézier segments and reduce each one (Eq. 5.41/5.42, end points fixed, pivot
   point averaged from both recursions),
2. join them C0 and restore the original continuity by knot removal (multiplicity s → s−1;
   simple knots stay simple),
3. **verify rigorously**: elevate the result back by one, refine both to common knots and compare
   control polygons. If the bound exceeds the tolerance with full continuity, the C0 version is
   checked; if that also fails, `GeometryError` is raised.

The bound is conservative, so a reduction may be rejected although the true deviation is
slightly below the tolerance. Undoing a degree elevation always succeeds exactly.

## Splitting and joining

- `split(u)` raises the multiplicity of u to p and cuts the control polygon; both parts are
  clamped and **keep their parameter intervals** ([a, u] and [u, b]).
- `join(a, b, tol)` unifies degrees (elevation), rescales b's weights so the joint weights agree,
  shifts b's interval to start at a's end, closes a gap ≤ tol at its midpoint and finally removes
  the joint knot as far as this is exact. Hence `join(*c.split(u))` reproduces `c.clamped()`.
- `join_u/join_v` additionally make the surfaces compatible across the seam (v degree, domain,
  knots). The seam rows must agree as **control rows**, i.e. both boundaries must share their
  parametrisation, not only their geometry; rational seams need proportional weights.

## Intersections

Curve–curve and curve–surface intersections (`analysis/intersect.hpp`) use subdivision of
Bézier pieces with bounding-box rejection. A pair of pieces becomes a leaf when it can contain
at most one intersection:

- curve/curve: the tangent cones (containing all control-polygon edges, hence all tangents)
  are separated, also against the reversed orientation – a second intersection would put the
  chord into both cones (it is a positive combination of tangents);
- curve/surface: the angle between curve tangents and surface normals stays below
  90° − (α + β), so the height of the curve over the surface is monotone.

Tangential contacts never separate the cones; there a size limit (1e-6 of the model size)
stops the subdivision. In each leaf, box-constrained Gauss–Newton on a(s) − b(t) = 0 or
c(t) − S(u, v) = 0 converges quadratically for zero residuals. Candidates are merged if they
are close or if the curve point at their parameter midpoint is still within tolerance (a
tangential contact keeps the gap below tol over a short stretch).

### Surface–surface intersection

`intersect(surface, surface)` (`analysis/surface_intersect.hpp`):

1. **Seeds.** Pairs of Bézier patches with overlapping boxes are subdivided. When their normal
   cones are separated, the pair contains no closed loop (Sederberg), so every intersection
   curve crosses the pair's boundary: the boundary iso-curves of each patch are intersected
   with the other patch (curve/surface intersection). Tangential regions stop at a size limit
   and seed by minimum-norm Gauss–Newton.
2. **Marching** along n1 × n2 (limit normals at poles): predictor in parameter space, corrector
   (S1 = S2 and the plane through the predicted point normal to the tangent), step control by
   the chord tolerance in 3D **and** in both parameter spaces (corrected vs. linear parameter
   midpoint). Branches end exactly on the domain boundary (one parameter fixed), close into
   loops (3D *and* parameters near the start – on closed surfaces the seam returns to the same
   3D point), or stop at tangential points.
3. **Knot lines.** Points are inserted exactly where a step crosses a knot line of either
   surface; curve and p-curves are interpolated piecewise between them (the parametrisation is
   only C1 there, e.g. the rational circle).
4. **Verified fit.** Between all fit points |S_i(pcurve_i(t)) − curve(t)| ≤ fit tolerance is
   checked and exact points are added where needed (near poles the p-curves vary quickly).
5. Seeds already on a traced curve are skipped – compared in 3D and in parameter distance
   weighted by |S_u|, |S_v| (seams: far apart; poles: u irrelevant).
6. **Tangential contacts.** A traced curve shorter than 1e-5 of the model size with parallel
   normals (marched around a contact point, e.g. on a seam) is reported as a point. A curve with
   nearly parallel normals everywhere (angle < 0.05) is tested at its middle: points of the
   first surface 1e-3 of the model size to either side of it lie on the same side of the second
   surface if the surfaces only touch (two cylinders side by side) – then it is reported as a
   point and further seeds on it are skipped – and on opposite sides if they cross.
7. **Touching the domain boundary.** A curve that touches a domain edge tangentially (a circle
   touching the edge of a plane's domain) ends where it comes within the tolerance of the edge
   and is continued from another seed. Open curves ending on the same boundary at coinciding
   points, or a short gap apart with tangents along the gap (sagitta ≤ 100 tol), are joined.

Curve portions lying in the other surface (e.g. an edge in a plane) are detected as overlaps
in the curve/curve and curve/surface intersections, which keeps such symmetric configurations
fast (plane through the middle of a torus along its knot lines: 0.4 s).

### Planar sections

`planar_section(surface, plane)` (`analysis/section.hpp`) intersects the surface with a
rectangle in the plane that covers the projection of the control points (plus margin), so no
section curve is missed; `pcurve2` and `u2, v2` are re-expressed in the plane's own orthonormal
(s, t) coordinates with origin at `plane.point`. `planar_sections` slices at equal spacing,
`intersect(curve, plane)` returns the crossing points of a curve.

## Trimmed surfaces

`TrimmedSurface` (`topology/trimmed_surface.hpp`) = surface + closed loops of curves in its
parameter domain (x = u, y = v, z = 0); the outer loop is oriented counter-clockwise, holes
clockwise. Helpers: `rectangle_loop`, `circle_loop` (exact), `polygon_loop`.

- `contains(u, v)`: even-odd test on dense polylines of the loops.
- `area()`: Green's theorem in the parameter plane - with F(u, v) = int |S_u x S_v| du the area
  is the sum of the loop integrals of F dv (holes subtract by orientation). Both integrals use
  the globally adaptive Gauss-Kronrod integrator with a bounded interval budget; the outer one
  has an absolute target, because along curves of constant v the integrand is rounding noise.
- `boundary_curves()`: trim curves mapped to 3D (cubic fit with verified tolerance, piecewise
  per trim-curve span).
- `split(surface, cutter)`: the p-curves of the surface/surface intersection cut the domain.
  Open curves (boundary to boundary - seams count as boundary) act as non-crossing chords of
  the parameter rectangle, giving #chords + 1 regions (traversal along the perimeter); closed
  curves become holes of the region around them and outer loops of new regions, also nested.

gmsh receives trimmed surfaces exactly: the loops are passed as B-spline wires in the
parameter space of the B-spline surface (`addBSplineSurface(..., wireTags, wire3D=False)`).

## Shells (B-rep)

`Shell::sew(faces, tol)` (`topology/shell.hpp`) maps the boundary loops of the trimmed faces to
3D and matches them: two boundary curves form the same edge if each lies on the other within
the sewing tolerance (7 samples, closest-point projection), in the same or opposite direction;
closed curves may start at different points. Degenerate boundary pieces (poles) are not edges;
the seam of a closed surface becomes an edge used twice by the same face. T-junctions
(partially overlapping curves) are resolved by splitting at the other curves' end points.

- `is_closed()`: every edge used exactly twice; `free_edges()`, `non_manifold_edges()`.
- `is_consistently_oriented()`: both uses of every shared edge run in opposite directions.
- `oriented()`: flips faces by breadth-first propagation over shared edges (throws for
  non-orientable configurations); closed shells end up with outward normals.
- `volume()`: divergence theorem, V = 1/3 Σ ∬ S·(S_u × S_v) du dv over the trimmed faces, with
  the same Green's-theorem integration as the trimmed area.

Trim curves are mapped to 3D piecewise between the knot lines of the surface they cross (the
rational circle is only C1 there; a single cubic per trim curve would need thousands of points).

## Bodies with several shells

A `Shell` may consist of several connected components (`Shell::components()`): the outer
boundary and the void boundaries of a hollow solid, or separate lumps (e.g. a difference that
cuts a body in two). `oriented()` first makes every component consistent (breadth-first over
shared edges) and then orients closed components by **nesting depth**: a component lying inside
an even number of other closed components (tested with `point_in_shell` from a point on it) is
an outer boundary with outward normals and positive volume; one inside an odd number bounds a
void, with normals pointing into the void and negative volume. The volume of the shell is then
the solid volume, Booleans accept hollow inputs, and the STEP writer exports voids as
`BREP_WITH_VOIDS`.

## Boolean operations

`boolean(A, B, op)` (`topology/boolean.hpp`) on closed shells:

1. Every pair of faces with overlapping boxes is intersected (surface/surface intersection);
   each intersection curve is clipped to the parameter intervals on which both p-curves lie in
   their trimmed faces (3D curve and p-curves share the parameter).
2. Every face is split along its cuts by a **planar arrangement** in its parameter domain
   (`split(TrimmedSurface, cuts)`): cuts clipped to the region, all curves split at their mutual
   intersections, dangling pieces removed, faces traced with "next half-edge clockwise" at every
   vertex; counter-clockwise cycles inside the face become regions, clockwise cycles become
   holes of the smallest enclosing region. Cuts may end on loops or other cuts and may cross.
3. Each piece is classified by `point_in_shell` (ray casting against the trimmed faces, rays
   that hit a boundary or graze a face are discarded, majority of three valid rays) at an
   interior point. Candidates come from inward rays from the outer loop, widest first; for a
   face that touches the other shell tangentially (the surface/surface intersection reported
   contact points) the first candidate at least 100 tolerances away from the touching faces is
   used – the widest one may lie on the contact line, where ray casting is undecided.
4. Selection – union: A out + B out; intersection: A in + B in; difference: A out + B in
   (flipped) – and sewing. The pieces keep the orientation of their outward-oriented shells, so
   inner voids stay inward (no re-orientation of the result).

Sewing splits boundary curves at the end points of other boundary curves lying on them
(T-junctions), e.g. a hole circle split at different points on the two faces sharing it.
**Coincident faces.** Before intersecting, each pair of faces is sampled (9×9 points of each
surface projected onto the other). Samples whose foot point lies on the other domain's boundary
are beyond its extent and ignored; if all remaining samples lie on the other surface, the faces
are *coincident* (if there are none, they only touch along their edges and are skipped). For
coincident faces the surface/surface intersection is degenerate and is not called: instead each
face is cut along the other face's loops, projected into its parameter domain (runs of the loop
within tolerance of the surface, end points by bisection, interpolated until the parameter error
is below 1e-8 of the domain). Cuts running along a loop or an earlier cut are dropped (the same
boundary often arrives twice, e.g. as a projected loop and as an intersection curve). A piece
whose interior point lies on a coincident face of the other shell is classified **on**, with the
same or the opposite outward normal:

| piece | union | intersection | difference A − B |
|---|---|---|---|
| A, on, same normal | keep | keep | drop |
| A, on, opposite normal | drop | drop | keep |
| B, on (either) | drop | drop | drop |

So two boxes sharing a face fuse into one solid, identical solids give themselves (union,
intersection) or nothing (difference), and an empty result is an empty shell.
**Tangential contacts** are not cut: touching faces are classified as a whole (see above), and
where an intersection curve touches a face boundary, short runs within 10 tolerances of the
boundary are taken as touching – inside runs are dropped, gaps between inside runs closed (the
in/out test flips arbitrarily there). Spheres touching at a point (outside, inside, on the
seam, at a pole), spheres and cylinders touching box faces, a cylinder inside a box touching
two sides and two cylinders side by side give exact volumes; the touching lumps of a union
stay separate components.

**Overlaps thinner than the tolerance** count as contacts: if the intersection curves of a face
pair bound a lens at most extent × angle / 4 thick (extent of the curves, largest angle between
the normals on them) and that is below the tolerance, the pair is not cut but treated as
touching (two spheres overlapping by 1e-7 give two touching lumps and an empty intersection).
Thin overlaps above the tolerance are cut: the arrangement finds an interior point of each
region by an inward ray (midpoint to the first crossing of the cycle), so slivers of 1e-5 –
the circle segments a box face cuts from a cylinder lid it overlaps by 1e-5 – are kept.

Not handled: coincident faces whose common part crosses the seam of a closed surface (the
projected loop would jump in the parameter domain).

## Fitting

| Algorithm | Function |
|---|---|
| Eq. 9.3–9.6 parameters | `data_parameters` (uniform, chord length, centripetal) |
| A9.1 global interpolation, knots Eq. 9.8 | `interpolate_curve` |
| Eq. 9.21/9.22 end derivatives | `interpolate_curve(..., end_derivatives)` (degree ≥ 2) |
| Eq. 9.63–9.69 least squares, ends interpolated | `approximate_curve`, `approximate_curve_tol` |
| A9.4 surface interpolation | `interpolate_surface` |
| row/column least squares (P&T 9.4.3) | `approximate_surface` |

Linear systems are solved with a small dense LU (partial pivoting) in `src/fitting/fit.cpp`,
factored once per parameter/knot set and reused for all rows/columns of a surface (ADR 0002:
no Eigen needed so far). `approximate_curve_tol` measures the residual at the data parameters
(an upper bound of the point-to-curve distance) and finds the smallest admissible number of
control points by exponential search and bisection.

## Loft (skinning)

`loft` follows P&T section 10.3: the sections are made compatible (`make_compatible`: clamped,
domain [0, 1], highest degree, union of knots), section parameters v_k are averaged over the
control-point columns (Eq. 10.8, centripetal/chord length) or given explicitly, and every column
is interpolated across the sections with degree min(q, K−1). Rational sections are interpolated
in homogeneous space – the interpolation is linear in the data, so w·P and w are interpolated
separately with identical knots. Linear data are reproduced exactly (equal circles → cylinder,
linearly growing circles → cone).

## Sweep

- `sweep_translational` is exact: tensor product with P_ij = P_i + Q_j − T(v0), w_ij = w_i w_j.
- `sweep` moves the profile (given at the trajectory start) with the **rotation-minimising
  frame**, computed by the double reflection method (Wang, Jüttler, Zheng, Liu 2008) on a table
  of 4096 steps plus one step to the query parameter; optional linear twist about the tangent.
- `sweep2` places the profile between two rails: translation with rail1, rotation so that the
  rail-to-rail direction is followed (orientation from rail1's tangent), uniform scaling by the
  rail distance.
- The swept surfaces are skinned (`loft`) through exactly placed profile copies. Because guide
  curves are only C1/C2 at their knots (e.g. the rational circle), the skin is built piecewise
  between the distinct knots of the guide curves – otherwise cubic skinning converges like h²
  instead of h⁴ – and the pieces are joined. Per piece the number of copies is doubled until the
  deviation from the ideal sweep at the midpoints is below `SweepOptions::tolerance`.

## Coons patches and Gordon surfaces

Both are Boolean sums S = L1 + L2 − T (P&T section 10.5): L1 skins the u-curves across v, L2
the v-curves across u, T interpolates the intersection points; all three are brought to common
degrees and knots and P = P1 + P2 − P3. The sum interpolates every network curve exactly, which
requires a **consistent network**: C_k(u_l) = D_l(v_k) with common intersection parameters
(determined by curve–curve intersection via alternating closest-point projection and checked).
A Coons patch is the 2 + 2 case (linear blending); its boundary orientations are fixed
automatically. The sum is only exact for non-rational curves, so rational curves are first
replaced by cubic interpolants at their original parameters, piecewise per knot span, within the
tolerance (quarter-cylinder test: radius error < 2·tol).

## Extrapolation

- **Natural:** the boundary span is made a Bézier piece (knot insertion), its homogeneous control
  points are re-computed for the longer interval by de Casteljau at t > 1, and the inserted knots
  are removed again where exact. The continuation is the same (rational) polynomial, hence C∞ at
  the old boundary, and conics/quadrics built from them stay exact (circle → circle, torus →
  torus). Non-positive extrapolated weights are rejected.
- **Linear:** a line along the end tangent whose parameter speed equals |C'(b)|, joined with
  `join` → C1 for non-rational geometry, G1 otherwise. For surfaces the iso-curves across the
  side must have constant weights.
- The extension length is an arc length (Newton on the length of the extension); surfaces measure
  it along the middle iso-curve. Original parameters are preserved.

## Offsets

The exact offset of a NURBS is not a NURBS in general (the unit normal contains a square root),
so offsets are fitted (`modeling/offset.hpp`). The parametrisation of the input is kept:
offset(t) approximates the exact offset point of c(t) (resp. S(u, v)), and the deviation from it is
the verified error.

- **Curves** (planar): o(t) = c(t) + d · normalize(n × c'(t)) for a curve in a plane with normal
  n (the control points must lie in that plane – equivalent to planarity of the curve, since the
  basis functions are linearly independent). With the unit tangent T and N = n × T,
  dN/ds = −κ_s T for the signed curvature κ_s = n · (c' × c'') / |c'|³, hence
  o'(t) = |c'|(1 − d κ_s) T: the offset is regular while d κ_s < 1 and reverses direction (cusp)
  beyond. `offset_curve` throws where d κ_s ≥ 1 (curvature from `analysis/properties.hpp`,
  checked at every fit and check point), where c' vanishes, and at tangent discontinuities that
  would open a gap > tol/2. Every Bézier segment of c is fitted separately with a cubic
  (least squares with given knots, p + 2 points per span, knots halved until 17 check points per
  span are within tol/2): rational curves such as circles are only C1 at their knots, and a cubic
  across such a joint converges like h² instead of h⁴. The pieces are joined with `join`
  (joint knots removed where exact).
- **Surfaces:** O(u, v) = S(u, v) + d · N(u, v) with the unit normal (limit normal at poles). With
  principal curvatures k_i relative to N, the offset's principal radii are 1/k_i − d; the offset
  degenerates where d k_i ≥ 1, which `offset_surface` rejects (principal curvatures on a 9 × 9
  grid per knot rectangle). Normal discontinuities across knot lines (creases) are rejected
  when they would open a gap > tol. The fit is a non-rational degree-5 tensor-product spline
  (error O(h⁶): about 10× fewer control points per direction than a bicubic at tol = 1e-7 for a
  sphere). The knot lines of S are break lines with multiplicity 6 − k for an input continuity
  C^k (offset C^(k−1)), clamped to [1, 5]. Knots are first inserted by 1D fits along iso-lines
  (5 per original span, tol/2 per direction), then the tensor fit (rows, then columns, as in
  `sag_surface`) is verified at 9 × 9 points per knot rectangle and failing rectangles are
  refined. Exact references in the tests: plane → plane, sphere → sphere R ± d, cylinder →
  cylinder, torus → tube radius r ± d (independently checked on 400 × 400 grids: the true error
  stays below tol).
- **Not handled:** global self-intersections (bottlenecks narrower than 2|d|) are neither
  detected nor removed, and the curvature check on a sample grid can miss a very narrow
  curvature peak between samples (the fit then usually fails to converge and throws).

## Fillets

Constant-radius fillets (`modeling/fillet.hpp`) use the offset construction: the centre of a
fillet of radius r has distance r from both inputs, so it lies on the intersection of their
offsets by ±r. The fitted offsets only provide start values; centres and contact points are then
solved exactly, so the fillet does not inherit the offset fitting error.

- **Curves** (`fillet_curves`, coplanar curves in a plane with normal n): the corner X is an
  intersection of a and b. Sign convention as for `offset_curve`: side +1 puts the centre to the
  left of the curve's direction of travel seen from the tip of n (offset by +r), −1 to the
  right. Automatic sides: on each curve the *kept part* is the side of X with the longer
  parameter interval (a polyline corner a → X → b keeps a before and b after X); the centre lies
  on the side of a towards b's kept part and vice versa. The start value (intersection of the
  two offset curves nearest to X) is refined by Gauss–Newton on
  a(s) + d_a N_a(s) = b(t) + d_b N_b(t) with N = normalize(n × c'), whose derivative is
  c'(1 − d κ_s) (see Offsets); s and t are then the exact foot points of the centre (tangent
  points). The arc is the exact rational arc (`make_arc`) from a(s) to b(t), oriented in the
  direction of travel of the trimmed a, so `trimmed_a`, `arc`, `trimmed_b` (the parts on the
  far side of the tangent points from X, reversed where necessary) form a G1 chain.
  References: two lines at an angle θ → setback r / tan(θ/2) from the corner, centre at
  r / sin(θ/2); line and circle → centre at height ±r over the line and distance R ± r from the
  circle centre (tests to 1e-10).
- **Surfaces** (`fillet_surfaces`, rolling ball): side_i = +1 puts the ball centre on the side
  of N_i = normalize(S_u × S_v) of surface i (offset by +r), −1 on the other side. The spine
  start values come from the surface–surface intersection of the two offset surfaces (the
  longest curve; offsets keep the parametrisation, so the intersection p-curves are also start
  values for the contact parameters on s1 and s2). At a spine parameter t the exact ball centre
  C is solved by Newton on |C − S1| = r, |C − S2| = r and (C − spine0(t)) · T = 0 (the plane
  through the start spine point normal to its tangent): the gradients of the distance functions
  are the unit vectors from the foot points (local closest-point projections, warm-started) to
  C. Since d|C(t) − S_i|/dt = 0, the true spine tangent is normal to both contact radii, so each
  section — the shorter circular arc from contact 1 to contact 2 around C — lies in the plane
  normal to the spine. Sections are exact rational arcs with a common number of ≤ 90° segments
  (A7.1, one knot vector for all), skinned with `loft` (cubic in v, homogeneous interpolation).
  The spine and the p-curves (u_i, v_i, 0) are cubic interpolations of the exact section data
  with the same parameters; the contact curves are the boundary iso-curves u = 0 and u = 1 of the
  fillet. Tolerance control: at t = 1/4, 1/2, 3/4 of every interval between sections, the
  fillet is compared with the exact section arc at 17 values of u, and spine, mapped p-curves
  S_i(pcurve_i(t)) with the exact centre and contacts; failing intervals are bisected (the mid
  check section becomes a new section) until all checks pass. References: two perpendicular
  planes → quarter cylinder (radius error < 1e-7, contact lines exact), plane and perpendicular
  cylinder, plane and sphere → |P − spine(t)| = r on all arcs, spine on the analytic circle,
  fillet normal parallel to the surface normals along the contact curves.
- **Not handled:** variable radius, fillets of more than one spine branch (the longest
  intersection curve is used), trimming of s1/s2 along the contact p-curves (the p-curves are
  returned for that), setbacks/vertex blends where several fillets meet, and curve fillets of
  curves that do not intersect (no extension).

## Integral and differential properties

- Arc length and area use adaptive Gauss–Kronrod 7/15 (QUADPACK nodes) on every knot span
  (tensor rule for areas); the integrands |C'| and |S_u × S_v| are smooth inside spans. Default
  relative accuracy `kIntegrationRelTol = 1e-10`.
- `parameter_at_length` inverts the arc length by safeguarded Newton (bisection fallback),
  accumulating the length incrementally between iterates.
- Surface curvature uses the fundamental forms with n = S_u × S_v / |S_u × S_v|; the sign of H,
  k1, k2 follows this orientation. At singular points (poles) the values are Richardson-
  extrapolated from two nearby regular points.

## Exactness notes

- Arcs use ≤ 90° segments with middle weight cos(Δθ/2); spheres, cylinders, cones and tori built by
  revolution are therefore exact up to floating point.
- Affine transforms act on control points only (affine invariance). Projective maps would also need
  weights and are not part of `Transform`.
- Poles (profile points on the axis) produce collapsed rows in the control net. There S_u (or S_v)
  vanishes and `normal()` returns the **limit normal** from the mixed derivative:
  S_u(u, v+h) ≈ h·S_uv, hence N ∝ ±S_uv × S_v, with the sign chosen so that it matches the
  neighbourhood (h points into the domain). If both partials vanish, `normal()` throws.

## Reference

L. Piegl, W. Tiller: *The NURBS Book*, 2nd ed., Springer, 1997.
