"""2D drawings of 3D parts: silhouettes, hidden lines (Make2D) and outlines."""
import math

import nurbforge as nf

# a plate with a through hole and a stud on top
plate = nf.boolean(nf.make_box((0, 0, 0), (4, 2, 1)),
                   nf.make_cylinder((1, 1, -1), (0, 0, 1), 0.5, 3), nf.BooleanOp.Difference)
stud = nf.make_cylinder((3, 1, 1), (0, 0, 1), 0.4, 1.5)
parts = [plate, stud]

# --- silhouettes: where a surface turns away from the viewer ---------------------------------
side = nf.silhouette(stud, direction=(0, 1, 0))  # parallel view along +y
print(f"silhouette of the stud along y: {len(side)} lines of length "
      f"{', '.join(f'{nf.length(c.curve):.3f}' for c in side)}")
ball = nf.make_sphere((0, 0, 0), 1.0)
rim = nf.silhouette(ball, eye=(3, 0, 0))  # perspective view from an eye point
print(f"sphere seen from distance 3: silhouette length {sum(nf.length(c.curve) for c in rim):.6f} "
      f"(exact {2 * math.pi * math.sqrt(1 - 1 / 9):.6f})")

# --- Make2D: visible and hidden lines in the drawing plane -----------------------------------
front = nf.make2d(parts, direction=(0, 1, 0), s_axis=(1, 0, 0))  # front view, x to the right
print(f"front view: {len(front.visible)} visible and {len(front.hidden)} hidden lines; "
      f"hidden length {sum(nf.length(c.curve) for c in front.hidden):.3f} (the hole)")
for c in front.hidden:
    a, b = c.curve(c.curve.domain[0]), c.curve(c.curve.domain[1])
    print(f"  hidden {c.kind.name.lower()} from ({a.x:.2f}, {a.y:.2f}) to ({b.x:.2f}, {b.y:.2f})")
iso = nf.make2d(parts, direction=(-1, -1.3, -0.8))
print(f"oblique view: {len(iso.visible)} visible, {len(iso.hidden)} hidden lines")

# --- visibility of any curves against the parts ----------------------------------------------
beam = nf.make_line((-1, 1, 0.5), (5, 1, 0.5))  # a line through the plate along x
pieces = nf.visibility([beam], parts, direction=(0, 0, -1))  # seen from above
print("line through the plate, seen from above: " +
      ", ".join(f"{'visible' if p.visible else 'hidden'} [{p.t0:.3f}, {p.t1:.3f}]"
                for p in pieces))

# --- outline: the contour of the parts seen along a direction --------------------------------
top = nf.outline(parts, direction=(0, 0, -1), s_axis=(1, 0, 0))
print(f"outline from above: {len(top.loops)} loops, area {nf.region_area(top.loops):.6f} "
      f"(exact {8 - math.pi * 0.25:.6f})")
