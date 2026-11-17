"""Hatching a planar region: lines, phase, minimum length, contours."""
import math

import nurbforge as nf

# region: 4 x 2 rectangle with a circular hole, as closed curve chains in the xy plane
# a closed polyline has kinks: split it into its straight pieces (needed for contours)
outer = nf.make_rectangle_curve((0, 0, 0), (4, 0, 0), (0, 2, 0)).decompose_bezier()  # CCW
hole = nf.make_circle((1.0, 1.0, 0), (0, 0, -1), 0.5)            # normal -z: clockwise
loops = [outer, [hole]]

segments, contours = nf.hatch(loops, (1, 0), 0.25)  # lines along x, 0.25 apart
print(f"{len(segments)} segments, total length {sum(s.length for s in segments):.6f}")
for s in segments[:6]:
    print(f"  line {s.line}: ({s.start.x:.4f}, {s.start.y:.4f}) -> ({s.end.x:.4f}, {s.end.y:.4f}), "
          f"loops {s.start_loop}->{s.end_loop}, region {s.region}")

# the lines are n . p = phase + k * spacing with n = (-d.y, d.x) = (0, 1): y = 0.1 + 0.25 k
shifted, _ = nf.hatch(loops, (1, 0), 0.25, phase=0.1)
print(f"phase 0.1: first line at y = {shifted[0].start.y:.4f}, {len(shifted)} segments")

# short pieces next to the hole can be dropped
kept, _ = nf.hatch(loops, (1, 0), 0.25, phase=0.1, min_length=0.6)
print(f"min_length 0.6: {len(kept)} of {len(shifted)} segments kept")

# diagonal hatching: direction (1, 1)
diag, _ = nf.hatch(loops, (1, 1), 0.25)
area = 8 - math.pi * 0.25
print(f"45 deg: {len(diag)} segments, spacing x length = {0.25 * sum(s.length for s in diag):.4f} "
      f"(region area {area:.4f})")

# contours: inward offsets by 0.1 and 0.3; at 0.3 the hole's contour merges with the outer one
_, cont = nf.hatch(loops, (1, 0), 0.25, contours=[0.1, 0.3])
for c in cont:
    L = sum(nf.length(k) for k in c.curves)
    print(f"contour d = {c.distance}: from loops {c.loops}, {len(c.curves)} pieces, length {L:.6f}")

# a neck narrower than 2 d splits a contour into two
dumbbell = nf.make_polyline([(0, 0, 0), (2, 0, 0), (2, 0.8, 0), (3, 0.8, 0), (3, 0, 0), (5, 0, 0),
                             (5, 2, 0), (3, 2, 0), (3, 1.2, 0), (2, 1.2, 0), (2, 2, 0), (0, 2, 0),
                             (0, 0, 0)])
for d in (0.1, 0.3):
    _, cs = nf.hatch([dumbbell.decompose_bezier()], (1, 0), 0.5, contours=[d])
    print(f"dumbbell, neck width 0.4, d = {d}: {len(cs)} contour(s)")
