"""Curve analysis: curvature, torsion, arc length, closest points and knot utilities."""
import math

import nurbforge as nf

# helix of radius 1 and pitch 2: curvature 1/(1 + c^2), torsion c/(1 + c^2) with c = pitch/(2 pi)
helix = nf.make_helix((0, 0, 0), (0, 0, 1), (1, 0, 0), pitch=2.0, turns=2.0, tolerance=1e-9)
c = 2.0 / (2 * math.pi)
print(f"curvature {nf.curvature(helix, 0.7):.8f}  exact {1 / (1 + c * c):.8f}")
print(f"torsion   {nf.torsion(helix, 0.7):.8f}  exact {c / (1 + c * c):.8f}")

# arc length of a circle and the inverse: the parameter at a given length
circle = nf.make_circle((0, 0, 0), (0, 0, 1), 2.0)
print(f"length {nf.length(circle):.12f}  (4 pi = {4 * math.pi:.12f})")
print(f"arc_length(0, 0.125) = {nf.arc_length(circle, 0.0, 0.125):.10f}")
u = nf.parameter_at_length(circle, math.pi)  # a quarter of the circumference
print(f"parameter at s = pi: u = {u:.12f}, point {circle(u)}")

# global closest point (point inversion) and a local search from a start value
cp = nf.closest_point(circle, (3, 3, 1))
print(f"closest point: u = {cp.u:.6f}, point {cp.point}, distance {cp.distance:.6f}")
print(f"expected distance sqrt((3 sqrt2 - 2)^2 + 1) = {math.hypot(3 * 2 ** 0.5 - 2, 1):.6f}")
for u0 in (0.8, 0.1):  # Newton from a start value; stops at the domain end, does not wrap
    loc = nf.closest_point_local(circle, (3, 3, 1), u0)
    print(f"local search from u0 = {u0}: u = {loc.u:.6f}, distance {loc.distance:.6f}")

# knot utilities
print("clamped_uniform_knots(3, 7):", nf.clamped_uniform_knots(3, 7))
print("unique_knots:", nf.unique_knots(circle.knots))
try:
    nf.validate_knots([0, 0, 1, 0.5, 1], 1, 3)
except nf.GeometryError as e:
    print("validate_knots:", e)

# make_compatible: common degree, knot vector and domain [0, 1]
line = nf.make_line((0, 0, 0), (1, 0, 0))
arc = nf.make_arc((0, 0, 0), (1, 0, 0), (0, 1, 0), 1.0, 0.0, math.pi)
for k in nf.make_compatible([line, arc]):
    print(f"compatible: degree {k.degree}, {len(k.control_points)} ctrl, domain {k.domain}")
