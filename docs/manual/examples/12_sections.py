"""Planar sections of solids: loops, plane coordinates, critical heights, extrusion bands."""
import math

import nurbforge as nf

# a plate with a bore and a spherical dome on top
plate = nf.make_box((-3, -2, 0), (6, 4, 1))
dome = nf.make_sphere((0, 0, 1), 1.5)
body = nf.boolean(plate, dome, nf.BooleanOp.Union)
body = nf.boolean(body, nf.make_cylinder((2, 0, -1), (0, 0, 1), 0.5, 4), nf.BooleanOp.Difference)
print(f"body: {len(body.faces)} faces, volume {body.volume(1e-8):.6f}")

# --- one section ------------------------------------------------------------
opts = nf.SectionOptions()
opts.s_axis = (1, 0, 0)  # first in-plane axis; t = normal x s
sec = nf.section(body, nf.Plane((0, 0, 0.5), (0, 0, 1)), opts)
print(f"origin {sec.origin}, s {sec.s_axis}, t {sec.t_axis}, normal {sec.normal}")
for i, loop in enumerate(sec.loops):
    kinds = ", ".join(f"face {c.face}{' exact' if c.exact else ''}" for c in loop.curves)
    print(f"loop {i}: area {loop.area:+.6f} ({'outer' if loop.area > 0 else 'hole'}), "
          f"{len(loop.curves)} piece(s): {kinds}")
print(f"material area {sum(l.area for l in sec.loops):.6f} (24 - pi/4 = {24 - math.pi / 4:.6f})")

# every piece has its curve in plane coordinates, in 3D and on its face, with the same parameter
c = sec.loops[1].curves[0]
t = 0.5 * sum(c.curve.domain)
p2, p3, n = c.curve(t), c.curve3d(t), c.normal(t)
print(f"hole piece at t: plane ({p2.x:.4f}, {p2.y:.4f}), 3D ({p3.x:.4f}, {p3.y:.4f}, {p3.z:.4f}), "
      f"outward normal ({n.x:+.3f}, {n.y:+.3f}, {n.z:+.3f})")

# --- through the dome -------------------------------------------------------
top = nf.section(body, nf.Plane((0, 0, 2), (0, 0, 1)))
r = math.sqrt(1.5**2 - 1.0)
print(f"z = 2: {len(top.loops)} loop, area {top.loops[0].area:.8f} (pi r^2 = {math.pi * r * r:.8f})")

# --- where does the topology of the sections change? ------------------------
hs = nf.critical_heights(body, (0, 0, 1))
print(f"critical heights: {[round(h, 6) for h in hs]}")
print(f"0.1..0.9 is an extrusion band: {nf.is_extrusion_band(body, (0, 0, 1), 0.1, 0.9)}; "
      f"1.1..2.0: {nf.is_extrusion_band(body, (0, 0, 1), 1.1, 2.0)}")
try:
    nf.section(body, nf.Plane((0, 0, 1.0), (0, 0, 1)))  # the plate top lies in this plane
except nf.GeometryError as exc:
    print(f"section at z = 1: GeometryError: {exc}")
