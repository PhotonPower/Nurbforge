"""Free-form surfaces through curves: loft, sweeps, Coons patch and Gordon surface."""
import math

import numpy as np
import nurbforge as nf

# --- loft through four circles (a vase); u along the sections, v across them
radii = [(0.0, 1.0), (1.0, 1.6), (2.5, 0.8), (3.5, 1.2)]  # (z, r)
sections = [nf.make_circle((0, 0, z), (0, 0, 1), r) for z, r in radii]
vase = nf.loft(sections, degree_v=3)
v_k = nf.loft_parameters(sections)
print(f"loft: net {vase.num_u} x {vase.num_v}, degree_v {vase.degree_v}, "
      f"section parameters {[round(v, 4) for v in v_k]}")
err = max(abs(math.hypot(vase(u, v).x, vase(u, v).y) - r)
          for (_, r), v in zip(radii, v_k) for u in np.linspace(0, 1, 9))
print(f"  surface passes through the sections: max radius error {err:.1e}")
uniform = nf.loft(sections, 3, [0.0, 1.0, 2.0, 3.0])  # explicit section parameters
print(f"  explicit v_params: domain_v {uniform.domain_v}, S(0, 2) = {uniform(0.0, 2.0)}")

# --- sweep a small circle along a quarter arc (rotation-minimising frame)
path = nf.make_arc((0, 0, 0), (1, 0, 0), (0, 1, 0), 5.0, 0.0, math.pi / 2)
start, tangent = path(0.0), path.tangent(0.0)
profile = nf.make_circle(tuple(start), tuple(tangent), 0.5)  # profile at the path start
bend = nf.sweep(profile, path, nf.SweepOptions(tolerance=1e-6))
err = max(abs(nf.closest_point(path, bend(u, v)).distance - 0.5)
          for u in np.linspace(0, 1, 9) for v in np.linspace(0.05, 0.95, 9))
print(f"sweep: net {bend.num_u} x {bend.num_v}, max deviation from tube radius {err:.1e}")
print(f"  area {nf.area(bend):.6f} (torus section 2 pi r * R pi/2 = "
      f"{2 * math.pi * 0.5 * 5 * math.pi / 2:.6f})")

# a twisted extrusion of a square: twist = total rotation at the end (radians)
square = nf.make_rectangle_curve((-0.5, -0.5, 0), (1, 0, 0), (0, 1, 0))
line = nf.make_line((0, 0, 0), (0, 0, 3))
twisted = nf.sweep(square, line, nf.SweepOptions(twist=math.pi / 2, tolerance=1e-5))
print(f"twisted sweep: corner at the top {twisted(0.0, 1.0)}")

# exact translational sweep: S(u, v) = C(u) + T(v) - T(v0)
wave = nf.interpolate_curve(np.array([(x, 0, 0.2 * math.sin(x)) for x in np.linspace(0, 6, 13)]))
slide = nf.sweep_translational(wave, nf.make_arc((0, 0, 0), (0, 1, 0), (0, 0, 1), 2.0, 0, 0.5))
print(f"translational sweep: degrees ({slide.degree_u}, {slide.degree_v}), "
      f"S(0, 1) - S(0, 0) = {slide(0, 1) - slide(0, 0)}")

# two-rail sweep: the profile is scaled so that it always spans the two rails
rail1 = nf.make_line((0, 0, 0), (0, 4, 0))
rail2 = nf.make_line((1, 0, 0), (2, 4, 0))
arch = nf.make_arc((0.5, 0, 0), (-1, 0, 0), (0, 0, 1), 0.5, 0.0, math.pi)  # from rail1 to rail2
fan = nf.sweep2(arch, rail1, rail2)
print(f"sweep2: span at the start {nf.distance(fan(0, 0), fan(1, 0)):.6f}, "
      f"at the end {nf.distance(fan(0, 1), fan(1, 1)):.6f}, top at the end {fan(0.5, 1).z:.6f}")

# --- Coons patch from four boundary curves forming a closed loop
c0 = nf.make_arc((1, 0, 0), (-1, 0, 0), (0, 0, 1), 1.0, 0.0, math.pi)  # S(u, 0)
c1 = nf.make_line((0, 2, 0), (2, 2, 0))                                # S(u, 1)
d0 = nf.make_line((0, 0, 0), (0, 2, 0))                                # S(0, v)
d1 = nf.make_line((2, 0, 0), (2, 2, 0))                                # S(1, v)
patch = nf.coons(c0, c1, d0, d1)
print(f"coons: S(0.5, 0) = {patch(0.5, 0)}, S(0.5, 0.5) = {patch(0.5, 0.5)}, "
      f"S(0.5, 1) = {patch(0.5, 1)}")

# --- Gordon surface through a consistent curve network
amps = {0.0: 0.2, 1.0: 0.5, 2.0: 0.3}  # amplitude of z = a sin(pi x / 2) at y
u_curves = [nf.interpolate_curve(np.array([(x, y, a * math.sin(math.pi * x / 2))
                                           for x in np.linspace(0, 2, 9)]),
                                 method=nf.Parametrization.Uniform) for y, a in amps.items()]
v_curves = [nf.interpolate_curve(np.array([(x, y, a * math.sin(math.pi * x / 2))
                                           for y, a in amps.items()]),
                                 degree=2, method=nf.Parametrization.Uniform)
            for x in (0.0, 1.0, 2.0)]
gordon = nf.gordon(u_curves, v_curves)
gap = max(nf.closest_point(gordon, c(t)).distance
          for c in u_curves + v_curves for t in np.linspace(0, 1, 11))
print(f"gordon: net {gordon.num_u} x {gordon.num_v}, network -> surface distance {gap:.1e}")
