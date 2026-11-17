"""Writing paths for two-photon lithography: slice a microlens, hatch every layer, export.

Units are micrometres. Every layer gets one contour (inward offset) plus parallel hatch lines
whose direction turns by 90 degrees from layer to layer; the hatch lines are written as a
meander (every second line reversed). The paths go to a plain text file, one polyline per
line: layer index, kind, then x y z triples.
"""
import math

import nurbforge as nf

# --- the part: plano-convex spherical cap, R = 30 um, sag 8 um, on a flat bas
R, H = 30.0, 8.0
a = math.sqrt(R * R - (R - H) ** 2)                  # aperture radius
centre = (0.0, 0.0, H - R)
meridian = nf.make_arc(centre, (0, 0, 1), (1, 0, 0), R, 0.0, math.asin(a / R))
cap = nf.revolve(meridian, centre, (0, 0, 1), 2 * math.pi)
base = nf.TrimmedSurface(nf.make_rectangle((-a, -a, 0), (2 * a, 0, 0), (0, 2 * a, 0)),
                         nf.circle_loop(0.5, 0.5, 0.5))
lens = nf.Shell.sew([nf.TrimmedSurface(cap), base]).oriented()
print(f"lens: aperture {2 * a:.3f} um, sag {H} um, volume {lens.volume():.4f} um^3, "
      f"critical heights {[round(h, 6) for h in nf.critical_heights(lens, (0, 0, 1))]}")

# --- process parameters -----------------------------------------------------
dz = 0.5          # layer distance (slicing)
spacing = 0.25    # hatch distance
inset = 0.15      # contour offset from the outer surface
opts = nf.SectionOptions()
opts.s_axis = (1, 0, 0)   # plane coordinates (s, t) = (x, y) in every layer


def polyline(curves, n=24):
    """Sample a chain of curves (plane coordinates) into one list of (s, t) points."""
    pts = []
    for c in curves:
        xy = c.sample(n)[:, :2]
        pts.extend(map(tuple, xy if not pts else xy[1:]))
    return pts


written, n_lines, n_contours = 0.0, 0, 0
with open("lens_paths.txt", "w") as out:
    out.write("# layer kind x y z ...  (micrometres)\n")
    for i in range(int(H / dz)):
        z = (i + 0.5) * dz                                  # mid-layer: never a critical height
        layer = nf.section(lens, nf.Plane((0, 0, z), (0, 0, 1)), opts)
        angle = (i % 2) * math.pi / 2
        segs, conts = nf.hatch(layer, (math.cos(angle), math.sin(angle)), spacing,
                               phase=0.5 * spacing, contours=[inset], min_length=0.05)
        for c in conts:
            pts = polyline(c.curves)
            out.write(f"{i} contour " + " ".join(f"{x:.4f} {y:.4f} {z:.4f}" for x, y in pts) + "\n")
        for k, s in enumerate(segs):                         # meander: reverse every second line
            p, q = (s.start, s.end) if k % 2 == 0 else (s.end, s.start)
            out.write(f"{i} hatch {p.x:.4f} {p.y:.4f} {z:.4f} {q.x:.4f} {q.y:.4f} {z:.4f}\n")
        length = sum(s.length for s in segs)
        written += spacing * dz * length
        n_lines += len(segs)
        n_contours += len(conts)
        if i % 5 == 0:
            print(f"layer {i:2d} z = {z:4.2f}: area {sum(l.area for l in layer.loops):8.3f} um^2, "
                  f"{len(segs):3d} hatch lines ({length:7.2f} um), {len(conts)} contour(s)")

cap_volume = math.pi * H * H * (3 * R - H) / 3
print(f"{int(H / dz)} layers, {n_lines} hatch lines, {n_contours} contours -> lens_paths.txt")
print(f"hatched volume {written:.2f} um^3 vs. cap volume {cap_volume:.2f} um^3")
