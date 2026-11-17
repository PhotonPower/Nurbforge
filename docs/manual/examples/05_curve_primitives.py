"""Curve primitives: lines, polylines, arcs, circles, polygons, rectangles and helices."""
import math

import nurbforge as nf

line = nf.make_line((0, 0, 0), (10, 0, 0))
poly = nf.make_polyline([(0, 0, 0), (4, 0, 0), (4, 3, 0), (0, 3, 0)])
arc = nf.make_arc((0, 0, 0), (1, 0, 0), (0, 1, 0), 5.0, 0.0, math.radians(120))
arc3 = nf.make_arc_3pt((5, 0, 0), (0, 5, 0), (-5, 0, 0))  # half circle through three points
circle = nf.make_circle(center=(0, 0, 0), normal=(0, 0, 1), radius=5.0)
circle3 = nf.make_circle_3pt((1, 0, 0), (0, 1, 0), (-1, 0, 0))
hexagon = nf.make_polygon((0, 0, 0), (0, 0, 1), 2.0, 6, x_axis=(1, 0, 0))
outline = nf.make_rectangle_curve((0, 0, 0), (4, 0, 0), (0, 2, 0))

for name, c in [("line", line), ("polyline", poly), ("arc 120 deg", arc), ("arc 3 pt", arc3),
                ("circle", circle), ("circle 3 pt", circle3), ("hexagon", hexagon),
                ("rectangle", outline)]:
    print(f"{name:12s} degree={c.degree} n_ctrl={len(c.control_points):2d} "
          f"rational={c.is_rational!s:5s} domain={c.domain} length={nf.length(c):.6f}")

# exactness: every point of the rational circle lies on the radius
r_err = max(abs(nf.norm(circle(i / 200)) - 5.0) for i in range(201))
print(f"circle radius error: {r_err:.1e}")

# helix with pitch 2 and 3 turns, a conical helix and a planar spiral (pitch 0)
helix = nf.make_helix((0, 0, 0), (0, 0, 1), (1, 0, 0), pitch=2.0, turns=3.0)
cone = nf.make_helix((0, 0, 0), (0, 0, 1), (1, 0, 0), pitch=2.0, turns=3.0, end_radius=0.25)
spiral = nf.make_helix((0, 0, 0), (0, 0, 1), (1, 0, 0), pitch=0.0, turns=4.0, end_radius=3.0)
print(f"helix: domain={helix.domain}, end point={helix(3.0)}, length={nf.length(helix):.6f}")
print(f"exact helix length 3*sqrt((2 pi)^2 + 2^2) = {3 * math.hypot(2 * math.pi, 2):.6f}")
print(f"conical helix end={cone(3.0)}, spiral end={spiral(4.0)}")
