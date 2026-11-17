"""Extending, offsetting and filleting surfaces."""
import math

import numpy as np
import nurbforge as nf

# --- extend: a cylinder band extended naturally stays a cylinder
band = nf.extrude(nf.make_arc((0, 0, 0), (1, 0, 0), (0, 1, 0), 1.0, 0.0, math.pi / 2), (0, 0, 1))
longer = nf.extend(band, nf.SurfaceSide.UMax, 0.5, nf.Extension.Natural)  # around the axis
taller = nf.extend(band, nf.SurfaceSide.VMax, 0.5, nf.Extension.Linear)   # along the axis
end = longer(longer.domain_u[1], 0.5)
print(f"extend UMax: domain_u {longer.domain_u}, end point {end}, radius {math.hypot(end.x, end.y):.12f}")
print(f"extend VMax: domain_v {taller.domain_v}, top z = {taller(0.5, taller.domain_v[1]).z:.6f}")

# --- offset_surface: sphere radius 1 offset by +0.5 along S_u x S_v
meridian = nf.make_arc((0, 0, 0), (0, 0, -1), (1, 0, 0), 1.0, 0.0, math.pi)
sphere = nf.revolve(meridian, (0, 0, 0), (0, 0, 1), 2 * math.pi)
n = sphere.normal(0.1, 0.5)
print(f"sphere normal at the equator {n} (points outward)")
off = nf.offset_surface(sphere, 0.5, tolerance=1e-6)
err = max(abs(nf.norm(off(u, v)) - 1.5) for u in np.linspace(0, 1, 21) for v in np.linspace(0, 1, 21))
print(f"offset +0.5: degrees ({off.degree_u}, {off.degree_v}), net {off.num_u} x {off.num_v}, "
      f"radius 1.5 within {err:.1e}")
try:
    nf.offset_surface(sphere, -1.2)  # inwards beyond the radius of curvature
except nf.GeometryError as e:
    print("offset -1.2 refused:", str(e)[:70], "...")

# --- fillet_surfaces: rolling ball between a floor and a wall
floor = nf.make_rectangle((-1, 0, 0), (5, 0, 0), (0, 3, 0))  # z = 0, normal +z
wall = nf.make_rectangle((0, 0, -1), (0, 3, 0), (0, 0, 5))   # x = 0, normal +x
sf = nf.fillet_surfaces(floor, wall, 0.5)  # side +1: ball on the +normal side of both
t = 0.5 * sum(sf.spine.domain)
print(f"fillet: net {sf.surface.num_u} x {sf.surface.num_v}, spine at mid {sf.spine(t)}")
print(f"  contact on floor {sf.contact1(t)}, on wall {sf.contact2(t)}")
q = sf.pcurve1(t)  # (u, v, 0) on the floor
print(f"  p-curve on the floor (u, v) = ({q.x:.6f}, {q.y:.6f}) -> {floor(q.x, q.y)}")
pts = sf.surface.sample(15, 15).reshape(-1, 3)
r_err = np.max(np.abs(np.hypot(pts[:, 0] - 0.5, pts[:, 2] - 0.5) - 0.5))
print(f"  fillet surface is the quarter cylinder r = 0.5 within {r_err:.1e}")
