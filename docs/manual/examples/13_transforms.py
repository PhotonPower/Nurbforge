"""Affine transforms of points, curves, surfaces and shells; patterns of copies."""
import math

import nurbforge as nf

T = nf.Transform.translation((10, 0, 0))
R = nf.Transform.rotation((0, 0, 0), (0, 0, 1), math.pi / 2)  # 90 deg about z (right hand)
S = nf.Transform.scaling((0, 0, 0), 2.0)                       # uniform
A = nf.Transform.scaling((0, 0, 0), (1, 2, 3))                 # per axis
M = nf.Transform.mirror((0, 0, 0), (1, 0, 0))                  # plane x = 0

p = nf.Vec3(1, 0, 0)
print("T p =", T.apply(p), " R p =", R.apply(p), " S p =", S.apply(p))
print("A (1,1,1) =", A.apply((1, 1, 1)), " M (2,3,4) =", M.apply((2, 3, 4)))
# composition: (T * R) applies R first, then T
print("(T * R) p =", (T * R).apply(p), " (R * T) p =", (R * T).apply(p))
# vectors are not translated
print("T.apply_vector((1,0,0)) =", T.apply_vector((1, 0, 0)))
print("matrix of T * R (3 x 4, row-major):", [round(x, 12) + 0.0 for x in (T * R).matrix()])
same = nf.Transform([0, -1, 0, 10, 1, 0, 0, 0, 0, 0, 1, 0])  # the same matrix by hand
print("hand-made matrix gives", same.apply(p))

# geometry is transformed exactly (affine invariance: only the control points move)
circle = nf.make_circle((0, 0, 0), (0, 0, 1), 1.0)
moved = circle.transformed(T * R)
print(f"circle moved: centre {(moved(0.0) + moved(0.5)) * 0.5}, length {nf.length(moved):.12f}")
ellipse = circle.transformed(A)
print(f"non-uniform scaling turns the circle into an ellipse: length {nf.length(ellipse):.8f}")

# shells: transformed(...) or the * operator; mirrors keep the normals outward
box = nf.make_box((0, 0, 0), (1, 2, 3))
for name, t in [("translated", T), ("rotated", R), ("scaled x2", S), ("mirrored", M)]:
    b = t * box
    print(f"box {name:10s}: volume {b.volume():8.4f}, oriented {b.is_consistently_oriented}")

# patterns: lists of transforms, applied with apply_pattern
hole = nf.make_circle((40, 0, 0), (0, 0, 1), 4.0)
ring = nf.apply_pattern(hole, nf.polar_pattern((0, 0, 0), (0, 0, 1), 8, 2 * math.pi))
print("polar 8 x 360 deg: centres at angles",
      [round(math.degrees(math.atan2(c(0.0).y + c(0.5).y, c(0.0).x + c(0.5).x))) for c in ring])
fan = nf.polar_pattern((0, 0, 0), (0, 0, 1), 3, math.pi / 2)
print("polar 3 x 90 deg: angles", [round(math.degrees(math.atan2(t.apply(p).y, t.apply(p).x)))
                                   for t in fan])
row = nf.linear_pattern((5, 0, 0), 4)
print("linear 4 x 5:", [t.apply((0, 0, 0)).x for t in row])
grid = nf.grid_pattern((5, 0, 0), 3, (0, 5, 0), 2)
print("grid 3 x 2:", [(t.apply((0, 0, 0)).x, t.apply((0, 0, 0)).y) for t in grid])
pillars = nf.apply_pattern(nf.make_cylinder((0, 0, 0), (0, 0, 1), 1.0, 4.0), grid)
print(f"{len(pillars)} pillars, total volume {sum(c.volume() for c in pillars):.6f} "
      f"(6 * 4 pi = {24 * math.pi:.6f})")
