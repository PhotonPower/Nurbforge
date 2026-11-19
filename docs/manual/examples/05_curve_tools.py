"""Curve tools: chamfers, tangent constructions, tween curves, best-fit line and plane,
projection onto a plane and curves on surfaces."""
import math

import nurbforge as nf

# --- chamfers: one corner, and every corner of a rectangle outline
a, b = nf.make_line((0, 0, 0), (4, 0, 0)), nf.make_line((4, 0, 0), (4, 3, 0))
ch = nf.chamfer_curves(a, b, 1.0, 0.5)
print(f"chamfer: {ch.point_a} -> {ch.point_b}, cut length {nf.length(ch.line):.6f}")
rect = nf.make_rectangle_curve((0, 0, 0), (4, 0, 0), (0, 2, 0))
cut = nf.chamfer_corners(rect, 0.5)
print(f"chamfered rectangle: length {nf.length(cut):.6f} "
      f"(12 - 8*0.5 + 4*0.5*sqrt 2 = {12 - 4 + 2 * math.sqrt(2):.6f})")

# --- tangent constructions: arc from start/tangent/end, tangents to circles
arc = nf.make_arc_tangent((0, 0, 0), (1, 0, 0), (-1, 1, 0))
print(f"arc_tangent: length {nf.length(arc):.6f} (3/4 of a unit circle = {1.5 * math.pi:.6f})")
pulley1 = nf.make_circle((0, 0, 0), (0, 0, 1), 2.0)
pulley2 = nf.make_circle((10, 0, 0), (0, 0, 1), 1.0)
belt = nf.tangent_lines(pulley1, pulley2)
print(f"common tangents of two circles: {len(belt)}, lengths "
      + ", ".join(f"{nf.length(t.line):.6f}" for t in belt))
from_point = nf.tangent_lines((6, 0, 0), pulley1)
print("tangents from (6, 0, 0): touching points "
      + ", ".join(f"({t.point_b.x:.6f}, {t.point_b.y:.6f})" for t in from_point))

# --- tween curves between a circle and a square outline
circle = nf.make_circle((0, 0, 0), (0, 0, 1), 1.0)
square = nf.make_polygon((0, 0, 0), (0, 0, 1), 2.0, 4)
tweens = nf.tween_curves(circle, square, 3)
print("tween lengths: " + ", ".join(f"{nf.length(c):.4f}" for c in tweens))

# --- best-fit line and plane through measured points
pts = [(x, y, 0.1 * x - 0.05 * y + 2 + 0.001 * (-1) ** (x + y)) for x in range(5) for y in range(5)]
pf = nf.fit_plane(pts)
n = pf.plane.normal
print(f"fit_plane: normal ({n.x:.6f}, {n.y:.6f}, {n.z:.6f}), max deviation "
      f"{pf.max_deviation:.2e}, rms {pf.rms:.2e}")
lf = nf.fit_line([(t, 2 * t + 0.01 * (-1) ** t, 1) for t in range(10)])
print(f"fit_line: direction ({lf.direction.x:.6f}, {lf.direction.y:.6f}, {lf.direction.z:.6f}),"
      f" max deviation {lf.max_deviation:.2e}")

# --- flatten a tilted circle onto the xy plane: an exact ellipse
tilted = nf.make_circle((0, 0, 5), (0, 1, 1), 2.0)
ellipse = nf.project_to_plane(tilted, nf.Plane((0, 0, 0), (0, 0, 1)))
ys = [ellipse(i / 100).y for i in range(101)]
print(f"projected circle: y range {max(ys) - min(ys):.6f} (2 * 2 cos 45 = {4 / math.sqrt(2):.6f})")

# --- curves on a cylinder: through points, and a geodesic offset of a ring
cyl = nf.extrude(nf.make_circle((0, 0, 0), (0, 0, 1), 1.0), (0, 0, 4))
pts = [(math.cos(t), math.sin(t), 0.5 + t) for t in (0.0, 0.5, 1.0, 1.5)]
on = nf.interpolate_on_surface(cyl, pts)
p = on.curve(0.5)
print(f"curve on the cylinder: |r| at t = 0.5 is {math.hypot(p.x, p.y):.9f}")
ring = nf.make_line((0, 0.25, 0), (1, 0.25, 0))  # p-curve of the circle at z = 1
up = nf.offset_on_surface(cyl, ring, 0.75)
print(f"geodesic offset of the ring by 0.75: z = {up.curve(0.3).z:.9f}")
