"""Joining, extending, offsetting and filleting curves."""
import math

import nurbforge as nf

# --- join: two curves that meet end to start become one curve
a = nf.make_line((0, 0, 0), (2, 0, 0))
b = nf.make_arc((2, 1, 0), (0, -1, 0), (1, 0, 0), 1.0, 0.0, math.pi / 2)  # tangent continuation
ab = nf.join(a, b)
print(f"join: degree {ab.degree}, domain {ab.domain}, length {nf.length(ab):.6f} "
      f"(2 + pi/2 = {2 + math.pi / 2:.6f})")

# --- extend: natural continuation keeps a circle a circle, linear adds a tangent line
arc = nf.make_arc((0, 0, 0), (1, 0, 0), (0, 1, 0), 1.0, 0.0, math.pi / 2)
nat = nf.extend(arc, nf.CurveEnd.End, 0.5, nf.Extension.Natural)
lin = nf.extend(arc, nf.CurveEnd.End, 0.5, nf.Extension.Linear)
print(f"extend natural: domain {nat.domain}, end {nat(nat.domain[1])}, "
      f"|end| = {nf.norm(nat(nat.domain[1])):.12f}")
print(f"extend linear : domain {lin.domain}, end {lin(lin.domain[1])}")
print(f"original parameters kept: {nf.distance(arc(0.5), nat(0.5)):.1e}")

# --- offset_curve: planar offset fitted to a tolerance, parameters kept
loop = nf.make_circle((0, 0, 0), (0, 0, 1), 1.0)
inner = nf.offset_curve(loop, 0.2, plane_normal=(0, 0, 1), tolerance=1e-7)  # left = inside
outer = nf.offset_curve(loop, -0.2, plane_normal=(0, 0, 1), tolerance=1e-7)
err = max(abs(nf.norm(inner(i / 400)) - 0.8) for i in range(401))
print(f"offset +0.2: radius 0.8 within {err:.1e}, degree {inner.degree}, "
      f"{len(inner.control_points)} control points")
print(f"offset -0.2: length {nf.length(outer):.8f} (2 pi 1.2 = {2 * math.pi * 1.2:.8f})")
try:
    nf.offset_curve(loop, 1.5)  # beyond the radius of curvature on the concave side
except nf.GeometryError as e:
    print("offset 1.5 refused:", str(e)[:60], "...")

# --- fillet_curves: exact arc between two intersecting curves
l1 = nf.make_line((0, 0, 0), (4, 0, 0))
l2 = nf.make_line((4, 0, 0), (4, 3, 0))
f = nf.fillet_curves(l1, l2, 0.5)
print(f"fillet: centre {f.center}, tangent points {f.point_a} / {f.point_b}")
print(f"        params on a/b: {f.param_a:.4f} / {f.param_b:.4f}, arc length "
      f"{nf.length(f.arc):.6f} (pi/4 = {math.pi / 4:.6f})")
chain = nf.join(nf.join(f.trimmed_a, f.arc, 1e-9), f.trimmed_b, 1e-9)
print(f"G1 chain length {nf.length(chain):.6f} (3.5 + pi/4 + 2.5 = {6 + math.pi / 4:.6f})")

# --- fillet_corners: round every corner of a polyline or polygon
rect = nf.make_rectangle_curve((0, 0, 0), (4, 0, 0), (0, 2, 0))
rounded = nf.fillet_corners(rect, 0.5)
print(f"rounded rectangle: length {nf.length(rounded):.6f} "
      f"(12 - 8*0.5 + 2 pi 0.5 = {12 - 4 + math.pi:.6f})")
