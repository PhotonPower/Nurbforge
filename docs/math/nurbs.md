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
| Tolerance | `nf::kTolerance = 1e-10` (model units) |

## Implemented algorithms (Piegl & Tiller numbering)

| Algorithm | Function |
|---|---|
| A2.1 FindSpan | `find_span` |
| A2.2 BasisFuns | `basis_functions` |
| A2.3 DersBasisFuns | `basis_function_derivatives` |
| A4.1 / A4.3 rational point | `NurbsCurve::evaluate`, `NurbsSurface::evaluate` |
| A4.2 RatCurveDerivs | `NurbsCurve::derivatives` |
| first-order RatSurfaceDerivs | `NurbsSurface::first_derivatives` |
| A5.1 CurveKnotIns | `NurbsCurve::insert_knot` |
| A7.1 MakeNurbsCircle | `make_arc`, `make_circle` |
| A8.1 MakeRevolvedSurf | `revolve` |

## Exactness notes

- Arcs use ≤ 90° segments with middle weight cos(Δθ/2); spheres, cylinders, cones and tori built by
  revolution are therefore exact up to floating point.
- Affine transforms act on control points only (affine invariance). Projective maps would also need
  weights and are not part of `Transform`.
- Poles (profile points on the axis) produce collapsed rows in the control net; normals there are
  undefined and `normal()` throws.

## Reference

L. Piegl, W. Tiller: *The NURBS Book*, 2nd ed., Springer, 1997.
