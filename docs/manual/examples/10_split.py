"""Splitting surfaces into trimmed regions: by a cutting surface or by parameter-space curves."""
import math

import nurbforge as nf

ball = nf.revolve(nf.make_arc((0, 0, 0), (0, 0, -1), (1, 0, 0), 2.0, 0.0, math.pi),
                  (0, 0, 0), (0, 0, 1), 2 * math.pi)

# --- split by a surface: a plane patch at z = 0.8 cut by the sphere R = 2 ---
plane = nf.make_rectangle((-3, -3, 0.8), (6, 0, 0), (0, 6, 0))
parts = nf.split(plane, ball)
for p in sorted(parts, key=lambda p: p.area(1e-8)):
    print(f"plane piece: {1 + len(p.holes)} loop(s), area {p.area(1e-8):.8f}")
print(f"disk pi (4 - 0.64) = {math.pi * (4 - 0.64):.8f}, frame 36 - disk = "
      f"{36 - math.pi * (4 - 0.64):.8f}")

# the other way round: the sphere cut by the plane -> cap and the rest of the sphere
pieces = nf.split(ball, plane)
caps = sorted(p.area(1e-8) for p in pieces)
h = 2 - 0.8
print(f"sphere pieces: {len(pieces)}, areas {caps[0]:.8f} + {caps[1]:.8f}; "
      f"cap 2 pi R h = {2 * math.pi * 2 * h:.8f}")

# no intersection: the surface comes back untrimmed
far = nf.make_rectangle((-3, -3, 5), (6, 0, 0), (0, 6, 0))
print(f"no intersection: {len(nf.split(far, ball))} piece, "
      f"untrimmed={nf.split(far, ball)[0].is_untrimmed}")

# --- split by parameter-space curves (x = u, y = v, z = 0) ------------------
square = nf.make_rectangle((0, 0, 0), (1, 0, 0), (0, 1, 0))
# open curves: chords from boundary to boundary that do not cross each other
chords = [nf.make_line((0.25, 0, 0), (0.25, 1, 0)),
          nf.make_line((0.6, 0, 0), (0.9, 1, 0))]
strips = nf.split(square, chords)
print(f"two chords: {len(strips)} regions, areas "
      f"{', '.join(f'{q.area():.4f}' for q in sorted(strips, key=lambda q: q.area()))}")
ring = nf.split(square, [nf.make_circle((0.5, 0.5, 0), (0, 0, 1), 0.3)])  # closed curve
print(f"closed curve: {len(ring)} regions, areas "
      f"{', '.join(f'{q.area():.6f}' for q in sorted(ring, key=lambda q: q.area()))} "
      f"(pi 0.09 = {math.pi * 0.09:.6f})")
