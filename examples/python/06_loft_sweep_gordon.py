"""Surfaces from curves: loft (skinning), sweep along a path, Coons and Gordon patches."""
import math

import numpy as np
import nurbforge as nf

# --- loft: a vase through four circular sections -------------------------------------------
radii = [(0.0, 1.0), (1.0, 1.6), (2.5, 0.8), (3.5, 1.2)]  # (z, r)
sections = [nf.make_circle((0, 0, z), (0, 0, 1), r) for z, r in radii]
vase = nf.loft(sections, degree_v=3)  # u along the sections, v across them
v = nf.loft_parameters(sections)      # v value of each section
err = max(abs(math.hypot(p.x, p.y) - r)
          for (_, r), vv in zip(radii, v) for p in (vase(u, vv) for u in np.linspace(0, 1, 9)))
print(f"loft: {vase.num_u} x {vase.num_v} control points, area {nf.area(vase):.4f}, "
      f"max radius error at the sections {err:.1e}")

# --- sweep: a circle along a helix with a rotation-minimising frame -------------------------
helix_pts = [(math.cos(a), math.sin(a), 0.3 * a) for a in np.linspace(0, 3 * math.pi, 40)]
helix = nf.interpolate_curve(np.array(helix_pts))
start, tangent = helix(0.0), helix.tangent(0.0)
profile = nf.make_circle(tuple(start), tuple(tangent), 0.15)  # profile at the path start
tube = nf.sweep(profile, helix, nf.SweepOptions(frame=nf.SweepFrame.RotationMinimizing,
                                                tolerance=1e-5))
# every point of the tube lies 0.15 from the helix
err = max(abs(nf.closest_point(helix, tube(u, w)).distance - 0.15)
          for u in np.linspace(0, 1, 7) for w in np.linspace(0.05, 0.95, 7))
print(f"sweep: {tube.num_u} x {tube.num_v} control points, "
      f"max deviation from tube radius 0.15: {err:.1e}")

# --- Coons patch from four boundary curves ----------------------------------------------------
# a half-circle arch at y = 0 and a lower, elliptic one at y = 2, joined by two straight edges
c0 = nf.make_arc((1, 0, 0), (-1, 0, 0), (0, 0, 1), 1.0, 0.0, math.pi)  # S(u, 0): x 0 -> 2
squash = nf.Transform.translation((0, 2, 0)) * nf.Transform.scaling((0, 0, 0), (1, 1, 0.5))
c1 = c0.transformed(squash)  # S(u, 1): half ellipse at y = 2
d0 = nf.make_line((0, 0, 0), (0, 2, 0))  # S(0, v)
d1 = nf.make_line((2, 0, 0), (2, 2, 0))  # S(1, v)
patch = nf.coons(c0, c1, d0, d1)
print(f"coons: S(0.5, 0) = {patch(0.5, 0.0)}, S(0.5, 0.5) = {patch(0.5, 0.5)}")

# --- Gordon surface through a curve network ---------------------------------------------------
# u-curves: z = a sin(pi x / 2) at y = 0, 1, 2 with amplitudes a; v-curves through the same
# points at x = 0, 1, 2 (the network must intersect consistently)
def u_curve(y, amp):
    xs = np.linspace(0, 2, 9)
    return nf.interpolate_curve(np.stack([xs, np.full_like(xs, y), amp * np.sin(math.pi * xs / 2)],
                                         axis=1), method=nf.Parametrization.Uniform)


amps = {0.0: 0.2, 1.0: 0.5, 2.0: 0.3}
u_curves = [u_curve(y, a) for y, a in amps.items()]
v_curves = []
for x in (0.0, 1.0, 2.0):
    pts = np.array([(x, y, a * math.sin(math.pi * x / 2)) for y, a in amps.items()])
    v_curves.append(nf.interpolate_curve(pts, degree=2, method=nf.Parametrization.Uniform))
gordon = nf.gordon(u_curves, v_curves)
gap = max(nf.closest_point(gordon, c(t)).distance
          for c in u_curves + v_curves for t in np.linspace(0, 1, 11))
print(f"gordon: {gordon.num_u} x {gordon.num_v} control points, "
      f"max distance network -> surface {gap:.1e}")
