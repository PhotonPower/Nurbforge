"""Exact surfaces from points and curves: bilinear patch, rectangle, extrusion, revolution,
ruled surface."""
import math

import nurbforge as nf

rect = nf.make_rectangle((0, 0, 0), (2, 0, 0), (0, 1, 0))  # S(u, v) = o + u*eu + v*ev
print(f"rectangle: degrees ({rect.degree_u}, {rect.degree_v}), net {rect.num_u} x {rect.num_v}, "
      f"S(0.5, 0.5) = {rect(0.5, 0.5)}, normal {rect.normal(0.5, 0.5)}, area {nf.area(rect):.10f}")

saddle = nf.make_bilinear((0, 0, 0), (1, 0, 0), (0, 1, 0), (1, 1, 1))  # p00, p10, p01, p11
print(f"bilinear: S(0.5, 0.5) = {saddle(0.5, 0.5)}, area {nf.area(saddle):.10f}")

# extrusion: u follows the profile, v the direction (degree 1)
tube = nf.extrude(nf.make_circle((0, 0, 0), (0, 0, 1), 1.0), (0, 0, 2))
print(f"extrude: degrees ({tube.degree_u}, {tube.degree_v}), net {tube.num_u} x {tube.num_v}, "
      f"area {nf.area(tube):.10f} (4 pi = {4 * math.pi:.10f})")

# revolution: u runs around the axis (rational quadratic), v follows the profile
meridian = nf.make_arc((0, 0, 0), (0, 0, -1), (1, 0, 0), 1.0, 0.0, math.pi)  # south -> north pole
sphere = nf.revolve(meridian, (0, 0, 0), (0, 0, 1), 2 * math.pi)
print(f"revolve: degrees ({sphere.degree_u}, {sphere.degree_v}), net "
      f"{sphere.num_u} x {sphere.num_v}, area {nf.area(sphere):.10f} (4 pi)")
print(f"  normal at the pole (limit normal): {sphere.normal(0.3, 0.0)}")
print(f"  normal at the equator: {sphere.normal(0.0, 0.5)}")

torus_profile = nf.make_circle((2, 0, 0), (0, 1, 0), 0.5)  # circle in the xz plane
torus = nf.revolve(torus_profile, (0, 0, 0), (0, 0, 1), 2 * math.pi)
print(f"torus area {nf.area(torus):.8f} (4 pi^2 R r = {4 * math.pi ** 2 * 2 * 0.5:.8f})")

# ruled surface between two curves: v = 0 on c0, v = 1 on c1
c0 = nf.make_circle((0, 0, 0), (0, 0, 1), 1.0)
c1 = nf.make_circle((0, 0, 1), (0, 0, 1), 0.5)
frustum = nf.ruled(c0, c1)
print(f"ruled: S(u, 1) radius {nf.norm(frustum(0.3, 1.0) - nf.Vec3(0, 0, 1)):.12f}, "
      f"area {nf.area(frustum):.8f} (pi (r0 + r1) s = "
      f"{math.pi * 1.5 * math.hypot(0.5, 1.0):.8f})")
