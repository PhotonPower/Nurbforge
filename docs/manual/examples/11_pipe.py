"""Solid pipes for flow and field simulations: a tapered elbow and a pipe with a bulge."""
import math

import nurbforge as nf

# an elbow: quarter circle of bend radius 3, pipe radius 0.4 tapering to 0.3
bend = nf.make_arc((0, 0, 0), (1, 0, 0), (0, 1, 0), 3.0, 0.0, math.pi / 2)
elbow = nf.make_pipe(bend, 0.4, end_radius=0.3)
print(f"elbow: closed={elbow.is_closed}, faces {sorted(elbow.face_labels)}")
for name in ("inlet", "outlet"):
    f = elbow.faces[elbow.faces_labeled(name)[0]]
    print(f"  {name}: area {f.area():.8f}")
print(f"  areas pi r^2: {math.pi * 0.16:.8f} and {math.pi * 0.09:.8f}")

# radius as a function of the normalised path parameter t in [0, 1]
line = nf.make_line((0, 0, 0), (4, 0, 0))
bulge = nf.make_pipe(line, lambda t: 0.5 + 0.2 * math.sin(math.pi * t))
exact = math.pi * 4 * (0.25 + 2 * 0.5 * 0.2 * 2 / math.pi + 0.04 / 2)
print(f"bulge: volume {bulge.volume():.8f} (pi int r^2 dx = {exact:.8f})")

# a closed path gives a ring (torus) without end faces
ring = nf.make_pipe(nf.make_circle((0, 0, 0), (0, 0, 1), 3.0), 0.4)
print(f"ring: {len(ring.faces)} face, volume {ring.volume():.8f} "
      f"(2 pi^2 R r^2 = {2 * math.pi**2 * 3 * 0.16:.8f})")
