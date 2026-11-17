"""Rotationally symmetric optical surfaces: exact conics and fitted even aspheres."""
import math

import numpy as np
import nurbforge as nf

# --- a sphere R = 40 (c = 1/40, k = 0) is represented exactly
sphere = nf.Asphere(curvature=1 / 40.0)
print(f"sphere: is_conic={sphere.is_conic}, max_radius={sphere.max_radius()}, "
      f"sag(10) = {sphere.sag(10.0):.12f} (exact {40 - math.sqrt(40**2 - 10**2):.12f})")
s_surf = nf.asphere_surface(sphere, 15.0)
r_err = max(abs(nf.norm(s_surf(u, v) - nf.Vec3(0, 0, 40)) - 40)
            for u in np.linspace(0, 1, 13) for v in np.linspace(0, 1, 13))
print(f"  surface: degrees ({s_surf.degree_u}, {s_surf.degree_v}), rational, "
      f"distance from the centre - R: {r_err:.1e}")

# --- conics: paraboloid (k = -1) and ellipsoid (k = -0.5)
for k in (-1.0, -0.5, 0.5):
    con = nf.Asphere(curvature=1 / 20.0, conic=k)
    surf = nf.asphere_surface(con, 10.0)
    err = max(abs(surf(u, v).z - con.sag(math.hypot(surf(u, v).x, surf(u, v).y)))
              for u in np.linspace(0, 1, 9) for v in np.linspace(0, 1, 9))
    print(f"conic k = {k:4.1f}: max_radius {con.max_radius():8.3f}, sag error {err:.1e}")

# --- an even asphere: conic plus A4 r^4 + A6 r^6, fitted to a tolerance
asph = nf.Asphere(curvature=1 / 40.0, conic=-0.8, even=[1e-6, -2e-9])
print(f"asphere: is_conic={asph.is_conic}, sag(15) = {asph.sag(15.0):.9f}")
for tol in (1e-6, 1e-9):
    prof = nf.asphere_profile(asph, 15.0, tolerance=tol)
    r = np.linspace(0, 15, 301)
    err = max(abs(nf.closest_point(prof, (x, 0, asph.sag(x))).distance) for x in r)
    print(f"  profile tol {tol:.0e}: degree {prof.degree}, {len(prof.control_points)} ctrl, "
          f"max distance to the sag {err:.1e}")
a_surf = nf.asphere_surface(asph, 15.0, tolerance=1e-8)
print(f"  surface domains u {a_surf.domain_u}, v {a_surf.domain_v} (v = r for fitted profiles)")
p = a_surf(0.1, a_surf.domain_v[1])  # the rim r = r_max
print(f"  rim point {p}, r = {math.hypot(p.x, p.y):.9f}, z - sag(15) = {p.z - asph.sag(15.0):.1e}")
print(f"  sphere surface domains u {s_surf.domain_u}, v {s_surf.domain_v} (exact conic)")

# odd terms: B1 r + B2 r^2 + ... (here a small cone term B1)
odd = nf.Asphere(curvature=0.0, odd=[0.01])
print(f"odd asphere sag(10) = {odd.sag(10.0):.6f}")
