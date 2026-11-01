"""Intersections and trimming: surface/surface intersection, holes, splitting, planar sections."""
import math

import nurbforge as nf


def sphere(radius):
    meridian = nf.make_arc((0, 0, 0), (0, 0, -1), (1, 0, 0), radius, 0.0, math.pi)
    return nf.revolve(meridian, (0, 0, 0), (0, 0, 1), 2 * math.pi)


ball = sphere(2.0)
cylinder = nf.extrude(nf.make_circle((0, 0, -3), (0, 0, 1), 1.0), (0, 0, 6))

# --- surface/surface intersection: sphere R=2 and cylinder r=1 -> two circles z = +-sqrt(3) ---
res = nf.intersect(ball, cylinder)  # SurfaceIntersection: .curves and isolated .points
# the circles come back in pieces, broken where they cross the seams of the closed surfaces
for sign in (-1, 1):
    pieces = [c for c in res.curves if sign * c.curve(0.5 * sum(c.curve.domain)).z > 0]
    L = sum(nf.length(c.curve) for c in pieces)
    print(f"intersection z = {sign * math.sqrt(3):+.6f}: {len(pieces)} pieces, total length "
          f"{L:.10f} (2 pi = {2 * math.pi:.10f})")
# each curve carries its p-curves (u, v, 0) on both surfaces
c = res.curves[0]
w = c.pcurve1(c.curve.domain[0])
print(f"p-curve on the sphere: (u, v) = ({w.x:.4f}, {w.y:.4f}) -> {ball(w.x, w.y)}")

# --- trimming: trim loops live in the (u, v) parameter space of the surface ----------------
square = nf.make_rectangle((0, 0, 0), (2, 0, 0), (0, 2, 0))  # u = x / 2, v = y / 2
holed = nf.TrimmedSurface(square).with_hole(nf.circle_loop(0.5, 0.5, 0.25))  # r = 0.5 in xy
print(f"plate with hole: area {holed.area():.10f} (4 - pi/4 = {4 - math.pi / 4:.10f}), "
      f"contains centre: {holed.contains(0.5, 0.5)}")
# on a 4 x 2 plate the same (u, v) circle is an ellipse with semi-axes 1 and 0.5
plate = nf.make_rectangle((0, 0, 0), (4, 0, 0), (0, 2, 0))
disk = nf.TrimmedSurface(plate, nf.circle_loop(0.5, 0.5, 0.25))  # outer loop only
print(f"elliptic disk on the 4 x 2 plate: area {disk.area():.10f} (pi/2 = {math.pi / 2:.10f})")

# --- split a plane along its intersection with the sphere -----------------------------------
plane = nf.make_rectangle((-3, -3, 0.8), (6, 0, 0), (0, 6, 0))
parts = nf.split(plane, ball)  # trimmed regions: the disk inside and the frame outside
areas = sorted(p.area(1e-8) for p in parts)
print(f"split: {len(parts)} parts, areas {areas[0]:.8f} + {areas[1]:.8f} = {sum(areas):.8f}; "
      f"disk pi (4 - 0.64) = {math.pi * (4 - 0.64):.8f}")

# --- planar sections: slice the sphere with five parallel planes -----------------------------
slices = nf.planar_sections(ball, nf.Plane((0, 0, -1.6), (0, 0, 1)), spacing=0.8, count=5)
for i, s in enumerate(slices):
    z = -1.6 + 0.8 * i
    L = sum(nf.length(cc.curve) for cc in s.curves)
    print(f"section z = {z:+.1f}: circumference {L:.6f} "
          f"(2 pi sqrt(4 - z^2) = {2 * math.pi * math.sqrt(4 - z * z):.6f})")
