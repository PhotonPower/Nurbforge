# Optical surfaces

`nurbforge/optics/surfaces.hpp` (Python: `nf.Asphere`, `nf.asphere_surface`, …) builds lens
surfaces as NURBS for meshing and downstream solvers.

## Conventions

- Optical axis **z**, vertex at the origin, sag z(x, y) along +z.
- Curvature c = 1/R (R > 0: centre of curvature on +z), conic constant k: k = 0 sphere,
  k = −1 paraboloid, k < −1 hyperboloid, −1 < k < 0 prolate and k > 0 oblate ellipsoid.
- Even asphere: z(r) = c r² / (1 + √(1 − (1+k) c² r²)) + A₄ r⁴ + A₆ r⁶ + …;
  odd terms B₁ r + B₂ r² + B₃ r³ + … are supported as well.
- Zernike polynomials with **Fringe** (unnormalised, Z4 = 2ρ² − 1, Z9 = 6ρ⁴ − 6ρ² + 1,
  Z16 = 6th-order spherical) or **Noll** indexing (normalised to unit RMS over the unit disk),
  ρ = r / norm_radius.
- XY polynomials: Σ c · (x/norm)^i (y/norm)^j. Zernike and XY sags are added to a base asphere.

## Representations and accuracy

| Surface | Representation | Accuracy |
|---|---|---|
| Sphere, ellipsoid, paraboloid, hyperboloid (`Asphere` without terms) | revolved **exact** conic: piecewise rational quadratic profile (≤ 45° tangent turn per segment, weight from the shoulder point) | exact up to rounding (implicit-equation residual ~1e-13), also hemispheres with vertical rim tangent |
| Even/odd asphere | cubic spline profile z(r) revolved exactly; x(t) = r stays exactly linear, so the error is a pure sag error; symmetric fit ⇒ slope 0 at the vertex | adaptive knot insertion until |Δz| ≤ tolerance at 8 check points per knot span |
| Zernike, XY polynomial, any `sag_surface(f, …)` | tensor-product spline graph over a rectangle; x, y exactly linear; polynomials of degree ≤ `degree` reproduced exactly | adaptive knot insertion until |Δz| ≤ tolerance at 8 × 8 check points per knot rectangle |

The tolerance is **verified at dense check points**, not proven between them; the tests
evaluate the surfaces at additional random points and find errors within 1.5–2× the tolerance.
For a stricter bound, pass a smaller tolerance.

Zernike and XY surfaces are graphs over a square (default `[-norm_radius, norm_radius]²`);
cutting them to a circular aperture needs trimming (roadmap 0.4). For steep conic bases the
square's corners may lie outside the conic's domain – `sag_surface` then reports the undefined
sag; reduce the half width.

## Example

```python
import nurbforge as nf
from nurbforge.mesh import mesh_surfaces

front = nf.asphere_surface(nf.Asphere(curvature=1/40, conic=-0.8, even=[1e-6, -2e-9]), 15.0)
back = nf.asphere_surface(nf.Asphere(curvature=-1/120), 15.0)
mesh_surfaces([front, back.transformed(nf.Transform.translation((0, 0, 6)))], "lens.msh", size=1.0)
```
