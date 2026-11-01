"""Freeform optics: asphere, Zernike and XY-polynomial sag surfaces checked against the formulas.

Optical axis z, vertex at the origin, curvature c = 1/R. See docs/optics.md for conventions.
"""
import math

import numpy as np
import nurbforge as nf

rng = np.random.default_rng(7)


def sag_error(surface, sag, half_width, n=200):
    """Max |z_surface - sag(x, y)| at random points; graph surfaces use u = x, v = y."""
    err = 0.0
    for x, y in rng.uniform(-half_width, half_width, size=(n, 2)):
        err = max(err, abs(surface(x, y).z - sag(x, y)))
    return err


# --- even asphere, fitted to 1e-8 (revolved spline profile) ---------------------------------
asph = nf.Asphere(curvature=1 / 40.0, conic=-0.8, even=[1e-6, -2e-9])
r_max = 15.0
a_surf = nf.asphere_surface(asph, r_max, tolerance=1e-8)
err = 0.0
for u, v in rng.uniform(0, 1, size=(200, 2)):  # revolved surface: compare z with sag(r)
    p = a_surf(u, v)
    err = max(err, abs(p.z - asph.sag(math.hypot(p.x, p.y))))
print(f"asphere: {a_surf.num_u} x {a_surf.num_v} control points, sag error {err:.1e} "
      f"(tolerance 1e-8), sag at rim {asph.sag(r_max):.6f}")

# analytic sag formula, written out
c, k, (a4, a6) = asph.curvature, asph.conic, asph.even
r = 10.0
z = c * r * r / (1 + math.sqrt(1 - (1 + k) * c * c * r * r)) + a4 * r**4 + a6 * r**6
print(f"asphere: sag(10) = {asph.sag(r):.12f}, formula {z:.12f}")

# a pure conic is exact: the vertex curvature of the surface equals c
sphere = nf.asphere_surface(nf.Asphere(curvature=1 / 50.0), 20.0)
k_vertex = nf.surface_curvature(sphere, 0.3, 0.0)  # v = 0 is the vertex (pole)
k_mid = nf.surface_curvature(sphere, 0.3, 0.6)
print(f"sphere R=50: principal curvatures at vertex ({abs(k_vertex.k1):.6f}, "
      f"{abs(k_vertex.k2):.6f}), at v=0.6 Gaussian {k_mid.gaussian:.6e} (1/R^2 = {1 / 2500:.6e})")

# --- Zernike freeform on a spherical base (Fringe indexing) --------------------------------
zern = nf.ZernikeSag(base=nf.Asphere(curvature=1 / 80.0), norm_radius=10.0,
                     indexing=nf.ZernikeIndexing.Fringe,
                     coeffs=[0, 0, 0, 0.002, 0.001, -0.0015, 0, 0, 0.0005])  # Z1..Z9
z_surf = nf.zernike_surface(zern, tolerance=1e-8)  # graph over [-10, 10]^2
print(f"zernike: {z_surf.num_u} x {z_surf.num_v} control points, "
      f"sag error {sag_error(z_surf, zern.sag, 10.0):.1e} (tolerance 1e-8)")
# Fringe Z4 = 2 rho^2 - 1 (defocus), Z9 = 6 rho^4 - 6 rho^2 + 1 (primary spherical)
rho = 0.5
print(f"zernike: Z4(0.5) = {nf.zernike(4, nf.ZernikeIndexing.Fringe, rho, 0.0):+.4f} "
      f"(2 rho^2 - 1 = {2 * rho**2 - 1:+.4f}), Noll j=11 -> (n, m) = "
      f"{nf.zernike_nm(11, nf.ZernikeIndexing.Noll)}")

# circular clear aperture of radius 9: trim the square graph (u = x, v = y)
aperture = nf.TrimmedSurface(z_surf, nf.circle_loop(0.0, 0.0, 9.0))
print(f"zernike: surface area of the trimmed aperture {aperture.area():.6f} "
      f"(slightly above the projected disk pi 9^2 = {math.pi * 81:.6f})")

# --- XY polynomial: c * (x/norm)^i * (y/norm)^j, here an astigmatic saddle + coma-like term ---
xy = nf.XYPolynomialSag(norm_radius=5.0, terms=[(2, 0, 0.05), (0, 2, -0.03), (1, 2, 0.01)])
xy_surf = nf.xy_polynomial_surface(xy, half_width=5.0)
x, y = 2.0, -1.5
formula = 0.05 * (x / 5) ** 2 - 0.03 * (y / 5) ** 2 + 0.01 * (x / 5) * (y / 5) ** 2
print(f"xy-poly: sag(2, -1.5) = {xy_surf(x, y).z:.15f}, formula {formula:.15f}, "
      f"max error {sag_error(xy_surf, xy.sag, 5.0):.1e} (polynomial reproduced exactly)")
k = nf.surface_curvature(xy_surf, 0.0, 0.0)
print(f"xy-poly: Gaussian curvature at the vertex {k.gaussian:.6e} "
      f"(saddle: 2*0.05/25 * 2*(-0.03)/25 = {(0.1 / 25) * (-0.06 / 25):.6e})")
