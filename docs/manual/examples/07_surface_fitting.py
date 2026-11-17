"""Fitting surfaces to a grid of points: interpolation and least squares."""
import numpy as np
import nurbforge as nf

# a 12 x 10 grid of z = sin(x) cos(y); shape (nu, nv, 3), first index along u
x, y = np.meshgrid(np.linspace(0, 2, 12), np.linspace(0, 1.5, 10), indexing="ij")
grid = np.stack([x, y, np.sin(x) * np.cos(y)], axis=2)
print("grid shape", grid.shape)

surf = nf.interpolate_surface(grid, degree_u=3, degree_v=3)
coarse = nf.approximate_surface(grid, degree_u=3, degree_v=3, ctrl_u=6, ctrl_v=5)

rng = np.random.default_rng(1)
samples = rng.uniform([0, 0], [2, 1.5], size=(50, 2))
for name, s in (("interpolate_surface", surf), ("approximate_surface", coarse)):
    dev = max(nf.closest_point(s, (a, b, np.sin(a) * np.cos(b))).distance for a, b in samples)
    corner = nf.distance(s(1.0, 1.0), nf.Vec3(*grid[-1, -1]))
    print(f"{name}: net {s.num_u} x {s.num_v}, max distance to z = sin(x) cos(y) {dev:.1e}, "
          f"corner error {corner:.1e}")

# the interpolant passes through every grid point (closest-point distance)
err = max(nf.closest_point(surf, tuple(grid[i, j])).distance for i in range(12) for j in range(10))
print(f"interpolant: max distance of the 120 grid points {err:.1e}")
