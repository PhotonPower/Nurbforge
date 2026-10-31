"""Basic curves and surfaces."""
import math

import nurbforge as nf

line = nf.make_line((0, 0, 0), (10, 0, 0))
circle = nf.make_circle(center=(0, 0, 0), normal=(0, 0, 1), radius=5.0)
arc = nf.make_arc((0, 0, 0), (1, 0, 0), (0, 1, 0), 5.0, 0.0, math.radians(120))

print(line, line(0.5))
print(circle, "rational:", circle.is_rational)
print("arc end:", arc(arc.domain[1]))
print("tangent of circle at u=0:", circle.tangent(0.0))
