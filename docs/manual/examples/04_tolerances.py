"""Tolerances: kernel defaults, per-operation Tolerance objects and verified fits."""
import math

import nurbforge as nf

t = nf.Tolerance()
print(f"defaults: length {t.length:g}, parametric {t.parametric:g}, angular {t.angular:g}")
print("repr rounds to 6 digits:", t)

loose = nf.Tolerance(length=1e-6)
c = nf.make_circle((0, 0, 0), (0, 0, 1), 1.0)
p = nf.closest_point(c, (3, 0, 0), loose)
print(f"closest point u = {p.u:.6f}, distance {p.distance:.12f}")

# an approximating construction: the result is checked against the request
helix = nf.make_helix((0, 0, 0), (0, 0, 1), (1, 0, 0), pitch=0.5, turns=3, tolerance=1e-8)
err = 0.0
for i in range(301):
    q = helix(i / 100)        # parameter = turns
    r = math.hypot(q.x, q.y)
    err = max(err, abs(r - 1.0), abs(q.z - 0.5 * i / 100))
print(f"helix: degree {helix.degree}, {len(helix.control_points)} control points, "
      f"max deviation {err:.1e} (tolerance 1e-8)")
