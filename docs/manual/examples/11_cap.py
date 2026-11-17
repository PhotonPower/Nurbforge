"""Closing the planar openings of a shell: an open pipe becomes a solid with annular end faces."""
import math

import nurbforge as nf

outer = nf.extrude(nf.make_circle((0, 0, 0), (0, 0, 1), 2.0), (0, 0, 3))
inner = nf.extrude(nf.make_circle((0, 0, 0), (0, 0, 1), 1.0), (0, 0, 3))
pipe = nf.Shell.sew([nf.TrimmedSurface(outer), nf.TrimmedSurface(inner).flipped()])
print(f"pipe wall: {len(pipe.faces)} faces, closed={pipe.is_closed}, "
      f"free edges {len(pipe.free_edges())}")

solid = nf.cap(pipe, label="end")
print(f"capped: {len(solid.faces)} faces, closed={solid.is_closed}, "
      f"volume {solid.volume():.8f} (pi (R^2 - r^2) h = {math.pi * 3 * 3:.8f})")
for i in solid.faces_labeled("end"):
    f = solid.faces[i]
    print(f"  end face: {len(f.holes)} hole, area {f.area():.8f} (3 pi = {3 * math.pi:.8f})")

# a non-planar opening cannot be capped
saddle = nf.Shell.sew([nf.TrimmedSurface(nf.make_bilinear((0, 0, 0), (1, 0, 1), (0, 1, 1),
                                                            (1, 1, 0)))])
try:
    nf.cap(saddle)
except nf.GeometryError as exc:
    print(f"saddle: GeometryError: {exc}")
