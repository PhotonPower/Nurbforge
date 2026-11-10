"""Slicing and hatching: writing paths of a plano-convex microlens for two-photon lithography.

The lens (spherical cap, units µm) is cut into layers with `section`; every layer is filled with
contour lines (inward offsets) and parallel hatch lines whose angle turns by 90° from layer to
layer. The written volume, spacing x layer distance x total hatch length, approaches the volume
of the cap as the spacings shrink.
"""
import math

import nurbforge as nf

# spherical cap: radius of curvature 30 µm, sag 8 µm on a flat base at z = 0
R, H = 30.0, 8.0
a = math.sqrt(R * R - (R - H) ** 2)  # aperture radius
centre = (0.0, 0.0, H - R)
rim = math.asin(a / R)  # polar angle of the rim
meridian = nf.make_arc(centre, (0, 0, 1), (1, 0, 0), R, 0.0, rim)  # pole -> rim in the xz plane
cap = nf.revolve(meridian, centre, (0, 0, 1), 2 * math.pi)
base = nf.TrimmedSurface(nf.make_rectangle((-a, -a, 0), (2 * a, 0, 0), (0, 2 * a, 0)),
                         nf.circle_loop(0.5, 0.5, 0.5))
lens = nf.Shell.sew([nf.TrimmedSurface(cap), base]).oriented()
print(f"lens: aperture {2 * a:.2f} µm, sag {H} µm, volume {abs(lens.volume()):.2f} µm^3")

dz, spacing, contour = 0.5, 0.2, 0.15  # layer distance, hatch distance, contour inset (µm)
written = 0.0
for i in range(int(H / dz)):
    z = (i + 0.5) * dz
    layer = nf.section(lens, nf.Plane((0, 0, z), (0, 0, 1)))
    angle = (i % 2) * math.pi / 2  # 0° and 90° in turn
    # contours vanish by themselves in layers narrower than 2 x the inset
    segments, contours = nf.hatch(layer, (math.cos(angle), math.sin(angle)), spacing,
                                  contours=[contour])
    length = sum(s.length for s in segments)
    written += spacing * dz * length
    if i % 4 == 0:
        print(f"z = {z:4.2f} µm: {len(segments):3d} hatch lines, {length:8.2f} µm, "
              f"{len(contours)} contour(s)")

cap_volume = math.pi * H * H * (3 * R - H) / 3
print(f"written volume {written:.1f} µm^3, cap volume {cap_volume:.1f} µm^3")
