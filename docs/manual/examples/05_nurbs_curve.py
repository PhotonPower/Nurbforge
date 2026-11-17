"""The NurbsCurve class: construction, evaluation and shape-preserving knot/degree operations."""
import numpy as np
import nurbforge as nf

# a cubic B-spline with 6 control points: knot vector of size n + p + 1 = 10
ctrl = np.array([[0, 0, 0], [1, 2, 0], [3, 3, 0], [5, 1, 0], [6, 3, 0], [8, 2, 0]], float)
knots = [0, 0, 0, 0, 0.3, 0.6, 1, 1, 1, 1]
c = nf.NurbsCurve(3, knots, ctrl)
print("degree", c.degree, "domain", c.domain, "clamped", c.is_clamped, "rational", c.is_rational)
print("C(0.5)  =", c(0.5))
d0, d1, d2 = c.derivatives(0.5, 2)  # point, first and second derivative
print("C'(0.5) =", d1, " C''(0.5) =", d2)
print("unit tangent", c.tangent(0.5))
pts = c.sample(5)  # (5, 3) array at equally spaced parameters
print("samples\n", np.round(pts, 4))

# the same control points with a uniform clamped knot vector
b = nf.NurbsCurve.from_points(3, ctrl)
print("from_points knots", [round(k, 4) for k in b.knots])

# knot insertion and refinement keep the shape
ins = c.insert_knot(0.45, 2)
ref = c.refine_knots([0.1, 0.2, 0.8, 0.9])
dev = max(nf.distance(c(t), ins(t)) + nf.distance(c(t), ref(t)) for t in np.linspace(0, 1, 101))
print(f"insert: {len(ins.control_points)} ctrl, refine: {len(ref.control_points)} ctrl, "
      f"max deviation {dev:.1e}")

# knot removal undoes the insertion exactly
back, removed = ins.try_remove_knot(0.45, 2)
print("removed", removed, "times ->", len(back.control_points), "ctrl")

# degree elevation is exact, degree reduction is verified against a tolerance
up = c.elevate_degree(2)
down = up.reduce_degree(1e-10).reduce_degree(1e-10)
print("elevated degree", up.degree, "n_ctrl", len(up.control_points),
      "-> reduced back to degree", down.degree)
try:
    c.reduce_degree(1e-6)  # a genuine cubic has no quadratic representation
except nf.GeometryError as e:
    print("reduce_degree refused:", str(e)[:70], "...")

# Bezier decomposition, splitting and sub-curves keep the parameter intervals
print("Bezier segments:", [s.domain for s in c.decompose_bezier()])
left, right = c.split(0.5)
print("split domains", left.domain, right.domain)
print("sub_curve domain", c.sub_curve(0.2, 0.7).domain)
r = c.reversed()
print("reversed: r(0) == c(1):", nf.distance(r(0.0), c(1.0)) < 1e-12)
print("reparametrized domain", nf.reparametrize(c, 0.0, 10.0).domain)

# a rational curve: the weights and the homogeneous control points
quarter = nf.NurbsCurve(2, [0, 0, 0, 1, 1, 1], np.array([[1, 0, 0], [1, 1, 0], [0, 1, 0]], float),
                        [1.0, 2 ** -0.5, 1.0])
print("quarter circle: weights", [round(w, 6) for w in quarter.weights],
      "|C(0.3)| =", round(nf.norm(quarter(0.3)), 12))
