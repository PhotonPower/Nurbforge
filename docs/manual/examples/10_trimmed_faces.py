"""Trim loops and trimmed surfaces: holes, areas, point classification, 3D boundaries."""
import math

import nurbforge as nf

# a 2 x 2 square: u = x / 2, v = y / 2
square = nf.make_rectangle((0, 0, 0), (2, 0, 0), (0, 2, 0))

# --- loops in the (u, v) plane ----------------------------------------------
rect = nf.rectangle_loop(0.1, 0.9, 0.1, 0.9)
circ = nf.circle_loop(0.5, 0.5, 0.25)
tri = nf.polygon_loop([(0.2, 0.2), (0.8, 0.2), (0.5, 0.7)])
for name, loop in (("rectangle", rect), ("circle", circ), ("triangle", tri)):
    print(f"{name:9s}: {len(loop.curves)} curve(s), signed area {loop.signed_area():+.6f}, "
          f"reversed {loop.reversed().signed_area():+.6f}")
print(f"circle loop: exact rational curve of degree {circ.curves[0].degree}, "
      f"{len(circ.polyline())} polyline points (16 per span)")

# --- trimmed faces ----------------------------------------------------------
plate = nf.TrimmedSurface(square)                       # untrimmed: outer loop = domain
holed = plate.with_hole(circ)                           # hole of radius 0.5 in xy
print(f"untrimmed: is_untrimmed={plate.is_untrimmed}, area {plate.area():.10f}")
print(f"with hole: area {holed.area():.10f} (4 - pi/4 = {4 - math.pi / 4:.10f})")
print(f"contains centre {holed.contains(0.5, 0.5)}, contains (0.1, 0.1) {holed.contains(0.1, 0.1)}")

# loops are oriented on construction: a clockwise outer loop is reversed automatically
frame = nf.TrimmedSurface(square, rect.reversed(), [circ])
print(f"outer signed area {frame.outer.signed_area():+.4f}, hole {frame.holes[0].signed_area():+.4f}, "
      f"area {frame.area():.10f} (4 * 0.64 - pi/4 = {4 * 0.64 - math.pi / 4:.10f})")

# --- boundaries in 3D -------------------------------------------------------
loops3d = holed.boundary_curves(1e-8)
for i, loop in enumerate(loops3d):
    L = sum(nf.length(c) for c in loop)
    print(f"3D loop {i} ({'outer' if i == 0 else 'hole'}): {len(loop)} curve(s), length {L:.8f}")
print(f"hole circumference 2 pi 0.5 = {math.pi:.8f}")

# --- orientation and labels -------------------------------------------------
n0 = holed.surface.normal(0.3, 0.3)
f = holed.flipped().with_label("cover")
n1 = f.surface.normal(0.3, 0.3)
print(f"normal {n0.z:+.0f} -> flipped {n1.z:+.0f}, label {f.label!r}, area kept {f.area():.10f}")

# --- trimming a curved surface: the same (u, v) circle on a cylinder --------
cyl = nf.extrude(nf.make_circle((0, 0, 0), (0, 0, 1), 1.0), (0, 0, 2))  # u around, v along z
window = nf.TrimmedSurface(cyl, nf.rectangle_loop(0.0, 0.25, 0.25, 0.75))
print(f"quarter of the cylinder, middle half: area {window.area():.10f} "
      f"(2 pi / 4 * 1 = {math.pi / 2:.10f})")
