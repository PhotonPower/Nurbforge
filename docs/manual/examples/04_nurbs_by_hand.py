"""A NURBS curve built from its definition: a quarter circle as a rational quadratic."""
import math

import nurbforge as nf

w = math.sqrt(0.5)                        # cos(45 deg): weight of the middle control point
quarter = nf.NurbsCurve(
    degree=2,
    knots=[0, 0, 0, 1, 1, 1],             # clamped: n_ctrl + degree + 1 = 6 knots
    control_points=[(1, 0, 0), (1, 1, 0), (0, 1, 0)],
    weights=[1.0, w, 1.0])
print(quarter, "rational:", quarter.is_rational, "clamped:", quarter.is_clamped)

# every point lies exactly on the unit circle
err = max(abs(nf.norm(quarter(i / 100)) - 1.0) for i in range(101))
print(f"max | |C(u)| - 1 | = {err:.2e}")

# the same curve without weights is a parabola, not a circle
parabola = nf.NurbsCurve(2, [0, 0, 0, 1, 1, 1], [(1, 0, 0), (1, 1, 0), (0, 1, 0)])
print(f"non-rational: |C(0.5)| = {nf.norm(parabola(0.5)):.6f}")

# knot insertion changes the representation, not the curve
refined = quarter.insert_knot(0.3, 2)
print("knots after inserting 0.3 twice:", refined.knots)
print(f"max deviation after insertion: "
      f"{max(nf.distance(quarter(i / 50), refined(i / 50)) for i in range(51)):.1e}")

# derivatives: C(u), C'(u), C''(u)
d = quarter.derivatives(0.0, 2)
print("C(0) =", d[0], " C'(0) =", d[1], " C''(0) =", d[2])
print(f"curvature at u = 0.5: {nf.curvature(quarter, 0.5):.12f}")

# invalid input raises GeometryError (a ValueError) with a precise message
try:
    nf.NurbsCurve(2, [0, 0, 1, 1], [(0, 0, 0), (1, 0, 0), (2, 0, 0)])
except nf.GeometryError as exc:
    print("GeometryError:", exc)
