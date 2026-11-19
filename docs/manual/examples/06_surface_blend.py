"""Blend and match surfaces: G0/G1/G2 transitions between surface sides."""
import math

import nurbforge as nf


def wavy(x0, y0, w, h, z0, phase):
    """Bicubic surface through z = f(x, y) on a 7 x 7 grid."""
    def z(x, y):
        return z0 + 0.3 * math.sin(1.3 * x + phase) * math.cos(0.9 * y)
    xs = [x0 + w * i / 6 for i in range(7)]
    ys = [y0 + h * j / 6 for j in range(7)]
    return nf.interpolate_surface([[(x, y, z(x, y)) for y in ys] for x in xs])


def join(b, v, s):
    """Largest gap, normal angle (degrees) and mean curvature jump along the iso-line v of b
    against the surface s."""
    gap = angle = jump = 0.0
    for k in range(1, 20):
        f = nf.closest_point(s, b(k / 20, v))
        cb, cs = nf.surface_curvature(b, k / 20, v), nf.surface_curvature(s, f.u, f.v)
        gap = max(gap, f.distance)
        angle = max(angle, math.degrees(math.asin(min(1.0, nf.norm(nf.cross(cb.normal,
                                                                             cs.normal))))))
        same = nf.dot(cb.normal, cs.normal) > 0
        jump = max(jump, abs(cb.mean - (cs.mean if same else -cs.mean)))
    return gap, angle, jump


a = wavy(0, 0, 2, 3, 0.0, 0.0)
b = wavy(3, -0.5, 2, 4, 0.8, 1.0)

# --- blend: x = 2 side of a (UMax) to x = 3 side of b (UMin)
for c in (0, 1, 2):
    bl = nf.blend_surfaces(a, nf.SurfaceSide.UMax, b, nf.SurfaceSide.UMin, continuity=c)
    (g0, a0, h0), (g1, a1, h1) = join(bl, 0.0, a), join(bl, 1.0, b)
    print(f"blend G{c}: degree across {bl.degree_v}, gap {max(g0, g1):.0e}, "
          f"normal angle {max(a0, a1):.0e} deg, curvature jump {max(h0, h1):.0e}")

# --- tension: how far the blend reaches along the surfaces before it turns
for t in (0.5, 1.0, 2.0):
    bl = nf.blend_surfaces(a, nf.SurfaceSide.UMax, b, nf.SurfaceSide.UMin, tension=t)
    print(f"tension {t}: middle cross-section length {nf.length(bl.iso_curve_u(0.5)):.4f}")

# --- exact blend from a cylinder rim to a flange (rational, constant weights across)
cyl = nf.revolve(nf.make_line((1, 0, 0), (1, 0, 2)), (0, 0, 0), (0, 0, 1), 2 * math.pi)
flange = nf.revolve(nf.make_line((2, 0, 3), (3, 0, 3)), (0, 0, 0), (0, 0, 1), 2 * math.pi)
neck = nf.blend_surfaces(cyl, nf.SurfaceSide.VMax, flange, nf.SurfaceSide.VMin)
radii = [math.hypot(neck(u / 8, 0.5).x, neck(u / 8, 0.5).y) for u in range(9)]
print(f"cylinder to flange: rational={bool(neck.weights)}, "
      f"radius at v = 0.5 varies by {max(radii) - min(radii):.1e}")

# --- match: move the VMin side of s onto the UMax side of a
s = wavy(2.3, -0.2, 2, 3.4, 0.4, 2.0).transposed()
for c in (0, 1, 2):
    m = nf.match_surface(s, nf.SurfaceSide.VMin, a, nf.SurfaceSide.UMax, continuity=c)
    g, ang, h = join(m, 0.0, a)
    far = max(nf.distance(m(u / 10, 1.0), s(u / 10, 1.0)) for u in range(11))
    print(f"match G{c}: gap {g:.0e}, normal angle {ang:.0e} deg, curvature jump {h:.0e}, "
          f"far side moved {far:.0e}")
