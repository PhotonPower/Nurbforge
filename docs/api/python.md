# Python API overview

```python
import nurbforge as nf
```

## Types

| Class | Key members |
|---|---|
| `Vec3(x, y, z)` | `x, y, z`, `+ - *`, iterable; tuples/lists convert implicitly |
| `NurbsCurve(degree, knots, control_points, weights=[])` | `evaluate(u)` / `c(u)`, `derivatives(u, order)`, `tangent(u)`, `sample(n) -> ndarray (n,3)`, `insert_knot(u, times)`, `refine_knots(X)`, `decompose_bezier()`, `reversed()`, `transformed(T)`, `domain`, `knots`, `weights`, `control_points`, `is_rational` |
| `NurbsCurve.from_points(degree, pts)` | clamped uniform B-spline |
| `NurbsSurface(du, dv, ku, kv, nu, nv, pts, weights=[])` | `evaluate(u,v)` / `s(u,v)`, `first_derivatives`, `normal`, `iso_curve_u/v`, `sample(nu,nv) -> ndarray (nu,nv,3)`, `refine_knots_u/v(X)`, `decompose_bezier()`, `transformed`, `transposed` |
| `Transform` | `translation`, `rotation(origin, axis, angle)`, `scaling`, `mirror(point, normal)`, `T1 * T2`, `apply`, `apply_vector` |

## Functions

| Group | Functions |
|---|---|
| Primitives | `make_line`, `make_polyline`, `make_arc`, `make_circle`, `make_bilinear`, `make_rectangle` |
| Construction | `extrude`, `revolve`, `ruled`, `reparametrize` |
| Arrays | `linear_pattern`, `grid_pattern`, `polar_pattern`, `apply_pattern` |
| Knots | `clamped_uniform_knots`, `unique_knots`, `validate_knots` |
| Vector math | `dot`, `cross`, `norm`, `distance` |

## Sub-packages

- `nurbforge.mesh` – `GmshSession`, `add_curve`, `add_surface`, `mesh_surfaces` ([meshing.md](../meshing.md))
- `nurbforge.llm` – `OllamaClient`, `GeometryAgent`, `Workspace`, CLI `nurbforge-llm` ([llm.md](../llm.md))

## Errors

`nurbforge.GeometryError` (subclass of `ValueError`) for invalid geometry.
