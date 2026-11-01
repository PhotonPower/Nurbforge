"""Fitting: interpolate and approximate measured points with curves and surfaces.

Interpolation passes exactly through the data; least-squares approximation trades exactness
for fewer control points. ``approximate_curve_tol`` picks the smallest control net that meets
a tolerance.
"""
import numpy as np
import nurbforge as nf

# --- curve data: 40 samples of a 3D spiral ------------------------------------------------
t = np.linspace(0.0, 1.0, 40)
pts = np.stack([np.cos(4 * t), np.sin(4 * t), 0.5 * t], axis=1)
params = nf.data_parameters(pts)  # centripetal parameters, the default of all fits


def max_dev(curve):
    """Largest distance between the data points and the curve at their parameters."""
    return max(nf.distance(curve(u), nf.Vec3(*p)) for u, p in zip(params, pts))


interp = nf.interpolate_curve(pts, degree=3)
print(f"interpolate_curve:     {len(interp.control_points):3d} control points, "
      f"max deviation {max_dev(interp):.1e}")

lsq = nf.approximate_curve(pts, degree=3, n_ctrl=8)
print(f"approximate_curve:     {len(lsq.control_points):3d} control points, "
      f"max deviation {max_dev(lsq):.1e}")

fit = nf.approximate_curve_tol(pts, degree=3, tolerance=1e-4)
print(f"approximate_curve_tol: {len(fit.control_points):3d} control points, "
      f"max deviation {max_dev(fit):.1e} (tolerance 1e-4)")

# prescribed end derivatives C'(0), C'(1) with respect to the curve parameter in [0, 1]:
# unit tangent of the spiral times its length (about the speed of a chord-length parameter)
L = float(np.sum(np.linalg.norm(np.diff(pts, axis=0), axis=1)))
d0 = L * np.array([0.0, 4.0, 0.5]) / np.hypot(4.0, 0.5)
d1 = L * np.array([-4 * np.sin(4.0), 4 * np.cos(4.0), 0.5]) / np.hypot(4.0, 0.5)
clamped = nf.interpolate_curve(pts, end_derivatives=(tuple(d0), tuple(d1)))
print(f"end derivatives prescribed: |C'(0) - d0| = "
      f"{np.linalg.norm(np.array(clamped.derivatives(0.0, 1)[1]) - d0):.1e}")

# --- surface data: 12 x 10 grid of z = sin(x) cos(y) ----------------------------------------
x, y = np.meshgrid(np.linspace(0, 2, 12), np.linspace(0, 1.5, 10), indexing="ij")
grid = np.stack([x, y, np.sin(x) * np.cos(y)], axis=2)  # shape (nu, nv, 3)

surf = nf.interpolate_surface(grid, degree_u=3, degree_v=3)
coarse = nf.approximate_surface(grid, degree_u=3, degree_v=3, ctrl_u=6, ctrl_v=5)

# deviation at random points of the underlying function: closest point on each surface
rng = np.random.default_rng(1)
for name, s in (("interpolate_surface", surf), ("approximate_surface", coarse)):
    dev = 0.0
    for xx, yy in rng.uniform([0, 0], [2, 1.5], size=(50, 2)):
        p = (xx, yy, np.sin(xx) * np.cos(yy))
        dev = max(dev, nf.closest_point(s, p).distance)
    print(f"{name}: {s.num_u} x {s.num_v} control points, "
          f"max distance to z = sin(x) cos(y): {dev:.1e}")
