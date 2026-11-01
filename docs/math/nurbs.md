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
