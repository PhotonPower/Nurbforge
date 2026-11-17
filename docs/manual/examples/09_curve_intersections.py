"""Curve/curve, curve/surface and curve/plane intersections."""
import math

import nurbforge as nf

# --- curve/curve: a line through a circle of radius 2 -> two crossing points 
circle = nf.make_circle((0, 0, 0), (0, 0, 1), 2.0)
line = nf.make_line((-3, 1, 0), (3, 1, 0))
for hit in nf.intersect(line, circle):  # sorted by the parameter s on the first curve
    print(f"line/circle: s = {hit.s:.6f}, t = {hit.t:.6f}, point = "
          f"({hit.point.x:+.6f}, {hit.point.y:+.6f}), gap {hit.distance:.1e}")
print(f"expected x = +-sqrt(3) = +-{math.sqrt(3):.6f}")

# a tangent line touches the circle once: the contact is found as long as the gap <= tol.length
tangent = nf.make_line((-3, 2, 0), (3, 2, 0))
print(f"tangent line: {len(nf.intersect(tangent, circle))} contact point(s)")

# skew lines do not intersect in 3D (they miss each other by 1)
skew = nf.make_line((0, -3, 1), (0, 3, 1))
print(f"skew lines: {len(nf.intersect(line, skew))} intersections")

# --- curve/surface: a helix piercing a plane patch several times ------------
helix = nf.make_helix((0, 0, 0), (0, 0, 1), (1, 0, 0), pitch=1.0, turns=3.0)
tilted = nf.make_rectangle((-2, -2, 0.0), (4, 0, 0), (0, 4, 3.0))  # z = 0.75 (y + 2)
hits = nf.intersect(helix, tilted)
print(f"helix/tilted plane: {len(hits)} points")
for h in hits:
    print(f"  t = {h.t:.6f} turns, (u, v) = ({h.u:.6f}, {h.v:.6f}), z = {h.point.z:.6f}")

# --- curve/plane: infinite plane, no patch needed ---------------------------
plane = nf.Plane((0, 0, 1.25), (0, 0, 1))
for h in nf.intersect(helix, plane):
    # u, v are coordinates along two orthonormal in-plane axes, origin at plane.point
    print(f"helix/plane z = 1.25: t = {h.t:.6f} (expected 1.25), point "
          f"({h.point.x:+.6f}, {h.point.y:+.6f}, {h.point.z:+.6f}), in-plane "
          f"radius {(h.u ** 2 + h.v ** 2) ** 0.5:.6f}")
print(f"signed distance of the origin from the plane: {plane.distance((0, 0, 0)):+.3f}")
