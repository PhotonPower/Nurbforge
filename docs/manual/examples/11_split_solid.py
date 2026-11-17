"""Splitting a solid by a plane: keep half of a symmetric part, with the cut faces labelled."""
import math

import nurbforge as nf

# a plate with a bore, symmetric about x = 0
part = nf.boolean(nf.make_box((-2, -1, 0), (4, 2, 1)),
                  nf.make_cylinder((0, 0, -1), (0, 0, 1), 0.4, 3), nf.BooleanOp.Difference)
left, right = nf.split(part, nf.Plane((0, 0, 0), (1, 0, 0)), cut_label="symmetry")
print(f"part volume {part.volume():.8f}, halves {left.volume():.8f} + {right.volume():.8f}")
print(f"right half: closed={right.is_closed}, labels {sorted(set(right.face_labels))}")
for f in (right.faces[i] for i in right.faces_labeled("symmetry")):
    n = f.surface.normal(0.5, 0.5)
    print(f"  symmetry face: area {f.area():.6f}, outward normal ({n.x:+.0f}, {n.y:+.0f}, {n.z:+.0f})")

# an oblique plane through a sphere: the two caps add up to the sphere
ball = nf.make_sphere((0, 0, 0), 1.0)
n = nf.Vec3(1, 2, 2) * (1 / 3)
below, above = nf.split(ball, nf.Plane(n * 0.3, n))
h = 0.7
print(f"cap above: {above.volume():.8f} (pi h^2 (3 - h) / 3 = {math.pi * h * h * (3 - h) / 3:.8f})")
print(f"rest below: {below.volume():.8f}, sum {below.volume() + above.volume():.8f} "
      f"(4 pi / 3 = {4 * math.pi / 3:.8f})")
