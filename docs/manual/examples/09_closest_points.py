"""Global and local closest-point queries (point projection / inversion)."""
import math

import nurbforge as nf

circle = nf.make_circle((0, 0, 0), (0, 0, 1), 2.0)
cp = nf.closest_point(circle, (3, 4, 0))
print(f"circle: u = {cp.u:.6f}, point ({cp.point.x:.6f}, {cp.point.y:.6f}), "
      f"distance {cp.distance:.6f} (expected 5 - 2 = 3)")

# local search: Newton iteration from a start parameter, finds the nearby (local) solution
line = nf.make_line((0, 0, 0), (10, 0, 0))
loc = nf.closest_point_local(line, (7.5, 1, 0), 0.5)
print(f"line, local from u0 = 0.5: u = {loc.u:.6f}, distance {loc.distance:.6f}")

ball = nf.revolve(nf.make_arc((0, 0, 0), (0, 0, -1), (1, 0, 0), 1.0, 0.0, math.pi),
                  (0, 0, 0), (0, 0, 1), 2 * math.pi)
sp = nf.closest_point(ball, (0, 3, 3))
print(f"sphere: (u, v) = ({sp.u:.6f}, {sp.v:.6f}), distance {sp.distance:.6f} "
      f"(expected 3 sqrt 2 - 1 = {3 * math.sqrt(2) - 1:.6f})")
# a point on the axis: every point of a parallel is closest -> one of them is returned
axis = nf.closest_point(ball, (0, 0, 5))
print(f"on the axis: distance {axis.distance:.6f}, point z = {axis.point.z:.6f}")
