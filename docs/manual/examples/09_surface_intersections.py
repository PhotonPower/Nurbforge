"""Surface/surface intersection, options, seams, planar sections."""
import math

import nurbforge as nf

ball = nf.revolve(nf.make_arc((0, 0, 0), (0, 0, -1), (1, 0, 0), 2.0, 0.0, math.pi),
                  (0, 0, 0), (0, 0, 1), 2 * math.pi)                    # sphere R = 2
tube = nf.extrude(nf.make_circle((0, 0, -3), (0, 0, 1), 1.0), (0, 0, 6))  # cylinder r = 1

res = nf.intersect(ball, tube)
print(f"sphere/cylinder: {len(res.curves)} curves, {len(res.points)} touching points")
for c in res.curves:
    a, b = c.curve.domain
    z = c.curve(0.5 * (a + b)).z
    # residual of the p-curves: S1(pcurve1(t)) and S2(pcurve2(t)) must lie on curve(t)
    err = 0.0
    for k in range(21):
        t = a + (b - a) * k / 20
        p1, p2 = c.pcurve1(t), c.pcurve2(t)
        err = max(err, nf.distance(ball(p1.x, p1.y), c.curve(t)),
                  nf.distance(tube(p2.x, p2.y), c.curve(t)))
    print(f"  z = {z:+.6f}, length {nf.length(c.curve):.6f}, closed={c.closed}, "
          f"{len(c.points)} traced points, p-curve residual {err:.1e}")
total = sum(nf.length(c.curve) for c in res.curves)
print(f"total length {total:.10f}, exact 2 * 2 pi = {4 * math.pi:.10f}")

# tighter options: smaller chord tolerance -> more traced points, smaller fit tolerance
opts = nf.SurfaceIntersectionOptions(chord_tolerance=1e-6, fit_tolerance=1e-10)
fine = nf.intersect(ball, tube, opts)
print(f"chord 1e-6: {sum(len(c.points) for c in fine.curves)} traced points "
      f"(default: {sum(len(c.points) for c in res.curves)})")

# two cylinders side by side touch along a line: reported as points, not as a curve
other = nf.extrude(nf.make_circle((2, 0, -3), (0, 0, 1), 1.0), (0, 0, 6))
touch = nf.intersect(tube, other)
print(f"touching cylinders: {len(touch.curves)} curves, {len(touch.points)} points")

# --- planar sections --------------------------------------------------------
cut = nf.planar_section(ball, nf.Plane((0, 0, 1), (0, 0, 1)))
c = cut.curves[0]
q = c.pcurve2(c.curve.domain[0])  # in-plane coordinates (s, t) about plane.point
print(f"section z = 1: {len(cut.curves)} curve(s), length {sum(nf.length(k.curve) for k in cut.curves):.8f} "
      f"(2 pi sqrt 3 = {2 * math.pi * math.sqrt(3):.8f}), first point in plane coords "
      f"radius {math.hypot(q.x, q.y):.8f}")
for i, s in enumerate(nf.planar_sections(ball, nf.Plane((0, 0, -1.5), (0, 0, 1)), 1.0, 4)):
    L = sum(nf.length(k.curve) for k in s.curves)
    z = -1.5 + i
    print(f"z = {z:+.1f}: {L:.8f} (exact {2 * math.pi * math.sqrt(4 - z * z):.8f})")
