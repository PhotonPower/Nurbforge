"""Thickening: a curved lens-holder shell and a washer get a wall thickness."""
import math

import nurbforge as nf

# a spherical dome (radius 10, opening angle 60 degrees) thickened by 0.5 along its normal
arc = nf.make_arc((0, 0, 0), (0, 0, 1), (1, 0, 0), 10.0, 0.0, math.radians(60))
dome = nf.Shell.sew([nf.TrimmedSurface(nf.revolve(arc, (0, 0, 0), (0, 0, 1), 2 * math.pi))])
n = dome.faces[0].surface.normal(0.5, 0.5)
p = dome.faces[0].surface(0.5, 0.5)
print(f"normal at ({p.x:.2f}, {p.y:.2f}, {p.z:.2f}): ({n.x:+.3f}, {n.y:+.3f}, {n.z:+.3f}) "
      f"- towards the centre, so the wall grows inwards")
shell = nf.thicken(dome, 0.5)
exact = 2 * math.pi / 3 * (10**3 - 9.5**3) * (1 - math.cos(math.radians(60)))
print(f"thick dome: closed={shell.is_closed}, faces {sorted(shell.face_labels)}, "
      f"volume {shell.volume():.6f} (spherical shell sector {exact:.6f})")

# a flat washer, symmetric about its mid-plane
sq = nf.make_rectangle((-2, -2, 0), (4, 0, 0), (0, 4, 0))
washer = nf.Shell.sew([nf.TrimmedSurface(sq, nf.circle_loop(0.5, 0.5, 0.5),
                                         [nf.circle_loop(0.5, 0.5, 0.25)])])
w = nf.thicken(washer, 0.4, symmetric=True)
print(f"washer: volume {w.volume():.8f} (pi (2^2 - 1^2) 0.4 = {math.pi * 3 * 0.4:.8f}), "
      f"walls {len(w.faces_labeled('wall'))}")

# creased shells cannot be thickened: the offsets of the faces do not meet
box = nf.make_box((0, 0, 0), (1, 1, 1))
tray = nf.Shell.sew([f for f in box.faces if f.label != "top"])
try:
    nf.thicken(tray, 0.1)
except nf.GeometryError as exc:
    print(f"open box: GeometryError: {exc}")
