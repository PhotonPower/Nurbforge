"""Points and vectors: Vec3, implicit conversions and vector helpers."""
import numpy as np
import nurbforge as nf

a = nf.Vec3(1.0, 2.0, 2.0)
b = nf.Vec3([0.0, 1.0, 0.0])          # from any sequence of three numbers
print("a + b =", a + b, " a - b =", a - b, " 2 * a =", 2 * a)
print("components:", a.x, a.y, a.z, " tuple:", a.to_tuple(), " list:", list(a))
print("dot", nf.dot(a, b), " cross", nf.cross(a, b), " norm", nf.norm(a),
      " distance", nf.distance(a, b))

# Wherever a point is expected, tuples, lists and NumPy vectors work as well
line1 = nf.make_line((0, 0, 0), [1, 0, 0])
pts = np.array([[0.0, 0.0, 0.0], [3.0, 4.0, 0.0]])
line2 = nf.make_line(pts[0], pts[1])  # rows of an (n, 3) array
print("length of line2:", nf.length(line2))

# Results come back as Vec3; convert explicitly for NumPy work
p = line2(0.5)
print("midpoint", p, "->", np.array(p.to_tuple()))

# Bulk evaluation returns NumPy arrays directly
samples = line2.sample(5)
print("sample(5):", samples.shape, samples[:, 0])
