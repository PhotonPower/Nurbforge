"""Offsets and fillets: offset curves and surfaces with verified tolerance, curve and surface
fillets."""
import math

import numpy as np
import nurbforge as nf

# --- offset curve: an exact ellipse (semi-axes 2 and 1), offset outwards by 0.2 -------------
loop = nf.make_circle((0, 0, 0), (0, 0, 1), 1.0).transformed(
    nf.Transform.scaling((0, 0, 0), (2, 1, 1)))  # counter-clockwise seen from +z
d = -0.2  # d > 0 offsets to the left of the curve direction seen from plane_normal (inwards)
off = nf.offset_curve(loop, d, plane_normal=(0, 0, 1), tolerance=1e-7)
# the offset is no longer an ellipse; check its distance to the ellipse at dense points
err = 0.0
for t in np.linspace(0, 1, 400):
    err = max(err, abs(nf.closest_point(loop, off(t)).distance - abs(d)))
print(f"offset_curve: degree {off.degree}, {len(off.control_points)} control points, "
      f"max |distance - 0.2| = {err:.1e} (tolerance 1e-7)")
try:
    nf.offset_curve(nf.make_circle((0, 0, 0), (0, 0, 1), 1.0), 1.5)  # inwards beyond the centre
except nf.GeometryError as exc:
    print(f"offset_curve(circle r=1, 1.5 inwards) -> GeometryError: {exc}")

# --- offset surface: a sphere R = 2 offset by 0.5 along its normal --------------------------
meridian = nf.make_arc((0, 0, 0), (0, 0, -1), (1, 0, 0), 2.0, 0.0, math.pi)
sphere = nf.revolve(meridian, (0, 0, 0), (0, 0, 1), 2 * math.pi)
sign = 1.0 if nf.dot(sphere.normal(0.3, 0.4), sphere(0.3, 0.4)) > 0 else -1.0  # outward?
shell = nf.offset_surface(sphere, sign * 0.5, tolerance=1e-6)
r = np.linalg.norm(shell.sample(25, 25), axis=-1)
print(f"offset_surface: degree ({shell.degree_u}, {shell.degree_v}), "
      f"{shell.num_u} x {shell.num_v} control points, radius {r.min():.7f}..{r.max():.7f} "
      f"(expected 2.5, tolerance 1e-6)")

# --- curve fillet: radius 0.5 between two lines meeting at a right angle -------------------
a = nf.make_line((-2, 0, 0), (0, 0, 0))
b = nf.make_line((0, 0, 0), (0, 2, 0))
f = nf.fillet_curves(a, b, 0.5)  # automatic side: the inner corner
print(f"fillet_curves: centre {f.center}, tangent points {f.point_a} / {f.point_b}")
chain = nf.join(nf.join(f.trimmed_a, f.arc, 1e-9), f.trimmed_b, 1e-9)  # G1 chain
print(f"fillet chain: length {nf.length(chain):.10f} "
      f"(1.5 + pi/4 + 1.5 = {3 + math.pi / 4:.10f}), arc rational: {f.arc.is_rational}")

# --- rolling-ball fillet between two perpendicular planes ----------------------------------
floor = nf.make_rectangle((-1, 0, 0), (5, 0, 0), (0, 3, 0))  # z = 0, normal +z
wall = nf.make_rectangle((0, 0, -1), (0, 3, 0), (0, 0, 5))   # x = 0, normal +x
sf = nf.fillet_surfaces(floor, wall, 0.5)  # ball on the +normal side of both: x > 0, z > 0
p = sf.surface.sample(15, 15)
dev = np.abs(np.hypot(p[..., 0] - 0.5, p[..., 2] - 0.5) - 0.5).max()
t = 0.5 * sum(sf.spine.domain)
q = sf.pcurve1(t)  # contact curve on the floor in its (u, v) parameters
print(f"fillet_surfaces: {sf.surface.num_u} x {sf.surface.num_v} control points, "
      f"max deviation from the r = 0.5 cylinder {dev:.1e}")
print(f"  spine at mid {sf.spine(t)}, contact on floor {sf.contact1(t)} = floor{(q.x, q.y)} "
      f"-> {floor(q.x, q.y)}")
