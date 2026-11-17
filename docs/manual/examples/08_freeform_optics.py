"""Freeform optical surfaces: Zernike and XY polynomial sags, arbitrary sag functions and
circular apertures."""
import math

import numpy as np
import nurbforge as nf

# --- Zernike polynomials and their indices
for j in (1, 4, 9, 11):
    n_f, m_f = nf.zernike_nm(j, nf.ZernikeIndexing.Fringe)
    n_n, m_n = nf.zernike_nm(j, nf.ZernikeIndexing.Noll)
    print(f"Z{j:<2d} Fringe (n, m) = ({n_f}, {m_f:2d}), Noll (n, m) = ({n_n}, {m_n:2d})")
print(f"Fringe Z4(rho=1) = {nf.zernike(4, nf.ZernikeIndexing.Fringe, 1.0, 0.0)}, "
      f"Noll Z4(rho=1) = {nf.zernike(4, nf.ZernikeIndexing.Noll, 1.0, 0.0):.6f} (sqrt 3)")

# --- a Zernike surface: base sphere R = 80 plus astigmatism, coma and spherical terms
zern = nf.ZernikeSag(base=nf.Asphere(curvature=1 / 80.0), norm_radius=10.0,
                     indexing=nf.ZernikeIndexing.Fringe,
                     coeffs=[0, 0, 0, 0.002, 0.001, -0.0015, 0, 0, 0.0005])  # Z1..Z9
z_surf = nf.zernike_surface(zern, tolerance=1e-8)  # graph over [-10, 10]^2, u = x, v = y
print(f"zernike surface: domain {z_surf.domain_u} x {z_surf.domain_v}, "
      f"net {z_surf.num_u} x {z_surf.num_v}, degree {z_surf.degree_u}")
rng = np.random.default_rng(0)
err = max(abs(z_surf(x, y).z - zern.sag(x, y)) for x, y in rng.uniform(-10, 10, (400, 2)))
print(f"  max |S_z - sag| at 400 random points: {err:.1e} (tolerance 1e-8)")

# --- an XY polynomial: terms (i, j, c) for c (x/norm)^i (y/norm)^j on a flat base
xy = nf.XYPolynomialSag(base=nf.Asphere(), norm_radius=5.0,
                        terms=[(2, 0, 0.01), (0, 2, -0.02), (2, 1, 0.003)])
xy_surf = nf.xy_polynomial_surface(xy, half_width=4.0, tolerance=1e-9)
print(f"xy polynomial: domain {xy_surf.domain_u}, net {xy_surf.num_u} x {xy_surf.num_v}, "
      f"S(2, 1).z = {xy_surf(2.0, 1.0).z:.12f}, sag(2, 1) = {xy.sag(2.0, 1.0):.12f}")

# --- any sag function z = f(x, y): here a sinusoidal grating on a rectangle
def grating(x, y):
    return 0.05 * math.sin(2 * math.pi * x / 2.5) * math.exp(-(y / 4) ** 2)


g_surf = nf.sag_surface(grating, -5, 5, -4, 4, degree=3, tolerance=1e-6)
err = max(abs(g_surf(x, y).z - grating(x, y)) for x, y in rng.uniform([-5, -4], [5, 4], (400, 2)))
print(f"sag_surface: net {g_surf.num_u} x {g_surf.num_v}, max error {err:.1e} (tolerance 1e-6)")

# --- circular clear aperture: trim the Zernike graph to radius 9 (u = x, v = y)
aperture = nf.TrimmedSurface(z_surf, nf.circle_loop(0.0, 0.0, 9.0))
print(f"aperture contains (8, 4): {aperture.contains(8.0, 4.0)}, (8, 5): {aperture.contains(8.0, 5.0)}")
print(f"aperture area {aperture.area():.6f} (flat disc pi 9^2 = {math.pi * 81:.6f})")

# --- curvature at the vertex of the base sphere
k = nf.surface_curvature(nf.zernike_surface(nf.ZernikeSag(base=nf.Asphere(curvature=1 / 80.0),
                                                          norm_radius=10.0)), 0.0, 0.0)
print(f"curvature at the vertex: k1 = {k.k1:.8f}, k2 = {k.k2:.8f} (1/80 = {1 / 80:.8f})")
