"""Fitting curves to measured points: interpolation, least squares, fit to a tolerance."""
import numpy as np
import nurbforge as nf

# 40 noisy-free samples of a 3D spiral, unevenly spaced (denser at the start)
t = np.linspace(0.0, 1.0, 40) ** 1.5
pts = np.stack([np.cos(4 * t), np.sin(4 * t), 0.5 * t], axis=1)

# the three parametrisations of the data points
for method in (nf.Parametrization.Uniform, nf.Parametrization.ChordLength,
               nf.Parametrization.Centripetal):
    u = nf.data_parameters(pts, method)
    print(f"{method.name:12s} parameters of points 0..3: {[round(x, 4) for x in u[:4]]}")


def max_dev(curve, method=nf.Parametrization.Centripetal):
    """Largest distance between the data points and the curve at their parameters."""
    u = nf.data_parameters(pts, method)
    return max(nf.distance(curve(ui), nf.Vec3(*p)) for ui, p in zip(u, pts))


def true_dev(curve):
    """Largest distance of the curve from the exact spiral (dense closest-point check)."""
    tt = np.linspace(0, 1, 400)
    return max(nf.closest_point(curve, (np.cos(4 * s), np.sin(4 * s), 0.5 * s)).distance
               for s in tt)


interp = nf.interpolate_curve(pts, degree=3)
print(f"interpolate_curve:      {len(interp.control_points):3d} ctrl, data deviation "
      f"{max_dev(interp):.1e}, distance to spiral {true_dev(interp):.1e}")

for n_ctrl in (6, 8, 12):
    lsq = nf.approximate_curve(pts, degree=3, n_ctrl=n_ctrl)
    print(f"approximate_curve n={n_ctrl:2d}: {len(lsq.control_points):3d} ctrl, data deviation "
          f"{max_dev(lsq):.1e}")

for tol in (1e-3, 1e-5):
    fit = nf.approximate_curve_tol(pts, degree=3, tolerance=tol)
    print(f"approximate_curve_tol {tol:.0e}: {len(fit.control_points):3d} ctrl, data deviation "
          f"{max_dev(fit):.1e}")

# end derivatives: C'(0) and C'(1) with respect to the curve parameter on [0, 1]
d0, d1 = (0.0, 4.0, 0.0), (0.0, 4.0, 0.0)
ends = nf.interpolate_curve(pts, degree=3, end_derivatives=(d0, d1))
print(f"with end derivatives: {len(ends.control_points)} ctrl (2 more), C'(0) = "
      f"{ends.derivatives(0.0, 1)[1]}")

# linear data are reproduced exactly by any degree
line_pts = np.array([[i, 2 * i, 0] for i in range(6)], float)
lin = nf.interpolate_curve(line_pts, degree=3, method=nf.Parametrization.ChordLength)
print(f"line data: C(0.5) = {lin(0.5)}")
