"""Defeaturing (filling small holes) and hollowing smooth solids."""
import math

import nurbforge as nf

# a plate with a bore (diameter 2), a small through hole and a small blind hole
part = nf.make_box((-3, -2, 0), (6, 4, 1))
part = nf.boolean(part, nf.make_cylinder((-1.5, 0, -1), (0, 0, 1), 0.2, 3), nf.BooleanOp.Difference)
part = nf.boolean(part, nf.make_cylinder((1.5, 0, 0.5), (0, 0, 1), 0.25, 2), nf.BooleanOp.Difference)
part = nf.boolean(part, nf.make_cylinder((0, 0, -1), (0, 0, 1), 1.0, 3), nf.BooleanOp.Difference)
clean = nf.remove_holes(part, max_size=0.8)  # holes smaller than 0.8 across
print(f"part: {len(part.faces)} faces, volume {part.volume():.8f}")
print(f"defeatured: {len(clean.faces)} faces, volume {clean.volume():.8f} "
      f"(24 - pi = {24 - math.pi:.8f})")

# coplanar pieces left by a union are merged into single faces
u = nf.boolean(nf.make_box((0, 0, 0), (2, 2, 1)), nf.make_box((1, 1, 0), (2, 2, 1)),
               nf.BooleanOp.Union)
m = nf.merge_faces(u)
top = m.faces[m.faces_labeled("top")[0]]
print(f"union of two boxes: {len(u.faces)} faces -> merged {len(m.faces)}, "
      f"top face area {top.area():.8f}, volume {m.volume():.8f}")

# a sphere hollowed to a wall of 0.1: an inner void
ball = nf.hollow(nf.make_sphere((0, 0, 0), 1.0), 0.1)
print(f"hollow sphere: {len(ball.components())} shells, volume {ball.volume():.8f} "
      f"(4/3 pi (1 - 0.9^3) = {4 / 3 * math.pi * (1 - 0.9**3):.8f})")

# a sphere made of two labelled hemispheres, opened at the top: a cup
upper = nf.revolve(nf.make_arc((0, 0, 0), (0, 0, 1), (1, 0, 0), 1.0, 0.0, math.pi / 2),
                   (0, 0, 0), (0, 0, 1), 2 * math.pi)
lower = nf.revolve(nf.make_arc((0, 0, 0), (1, 0, 0), (0, 0, -1), 1.0, 0.0, math.pi / 2),
                   (0, 0, 0), (0, 0, 1), 2 * math.pi)
sphere = nf.Shell.sew([nf.TrimmedSurface(upper).with_label("top"),
                       nf.TrimmedSurface(lower).with_label("bottom")]).oriented()
cup = nf.hollow(sphere, 0.1, open_labels=["top"])
print(f"cup: closed={cup.is_closed}, labels {sorted(set(cup.face_labels))}, "
      f"volume {cup.volume():.8f} (2/3 pi (1 - 0.9^3) = {2 / 3 * math.pi * (1 - 0.9**3):.8f})")
