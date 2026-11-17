"""2D profiles: planar region Booleans, blend curves between curve ends, projection onto a surface."""
import math

import nurbforge as nf

# --- a flange profile: a plate with rounded ends minus four bolt holes and a central bore
plate = nf.make_rectangle_curve((-4, -1.5, 0), (8, 0, 0), (0, 3, 0))
ends = [[nf.make_circle((x, 0, 0), (0, 0, 1), 1.5)] for x in (-4, 4)]
profile = [[plate]]
for e in ends:
    profile = nf.region_boolean(profile, [e], nf.BooleanOp.Union)
holes = [[nf.make_circle((x, 0, 0), (0, 0, 1), 0.5)] for x in (-4, 4)]
holes.append([nf.make_circle((0, 0, 0), (0, 0, 1), 1.0)])
for h in holes:
    profile = nf.region_boolean(profile, [h], nf.BooleanOp.Difference)
exact = 8 * 3 + math.pi * 1.5**2 - 2 * math.pi * 0.25 - math.pi
print(f"flange: {len(profile)} loops, area {nf.region_area(profile):.8f} (exact {exact:.8f})")
outer = max(profile, key=lambda loop: nf.region_area([loop]))
print(f"outer loop: {len(outer)} exact curves (lines and arcs)")

# --- a G2 blend between the end of a straight lead-in and the start of an arc
lead = nf.make_line((-3, 2, 0), (-1, 2, 0))
arc = nf.make_arc((2, 0, 0), (1, 0, 0), (0, 1, 0), 1.0, math.pi / 2, math.pi)
blend = nf.blend_curves(lead, nf.CurveEnd.End, arc, nf.CurveEnd.Start, continuity=2)
print(f"blend: degree {blend.degree}, curvature at its ends "
      f"{nf.curvature(blend, 0.0):.6f} (line: 0) and {nf.curvature(blend, 1.0):.6f} (arc: 1)")

# --- an aperture outline projected onto a freeform lens surface along the optical axis
#     (a graph surface z(x, y) with u = x, v = y; revolved surfaces have a seam - split curves
#     that cross it)
lens = nf.zernike_surface(nf.ZernikeSag(base=nf.Asphere(curvature=1 / 20), norm_radius=6.0,
                                        coeffs=[0, 0, 0, 0.01, 0.005]))
outline = nf.make_circle((0.5, 0, 10), (0, 0, 1), 3.0)  # off-centre circle above the lens
on_lens, uv = nf.project_curve(outline, lens, direction=(0, 0, -1), tolerance=1e-7)
p = on_lens(0.25)
print(f"projected point ({p.x:.6f}, {p.y:.6f}, {p.z:.6f}); "
      f"sag there {lens(uv(0.25).x, uv(0.25).y).z:.6f}")
