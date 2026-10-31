"""Array functions: polar array of holes (circles) on a flange."""
import math

import nurbforge as nf

hole = nf.make_circle((40, 0, 0), (0, 0, 1), 4.0)
holes = nf.apply_pattern(hole, nf.polar_pattern((0, 0, 0), (0, 0, 1), 8, 2 * math.pi))
for i, h in enumerate(holes):
    c = h(0.0)
    print(f"hole {i}: start point ({c.x:7.3f}, {c.y:7.3f})")

grid = nf.apply_pattern(nf.make_rectangle((0, 0, 0), (1, 0, 0), (0, 1, 0)),
                        nf.grid_pattern((1.5, 0, 0), 3, (0, 1.5, 0), 2))
print(len(grid), "patches in grid")
