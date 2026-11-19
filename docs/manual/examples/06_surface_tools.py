"""Surface tools: sweeps with several profiles, rail revolve, tapered extrusion, ribbon, cone to a
point, exact Coons patches from arcs, three-sided patches and rebuilding."""
import math

import nurbforge as nf

# --- sweep with several profiles: a duct from a circle to a larger circle and back
path = nf.make_line((0, 0, 0), (0, 0, 6))
profiles = [nf.make_circle((0, 0, z), (0, 0, 1), r) for z, r in ((0, 1.0), (3, 1.5), (6, 1.0))]
duct = nf.sweep(profiles, path)  # parameters: where the profiles cross the path
mid = duct(0.25, 0.5)
print(f"duct: radius at z = 3: {math.hypot(mid.x, mid.y):.6f}, at the ends "
      f"{math.hypot(duct(0.25, 0).x, duct(0.25, 0).y):.6f}")

# --- rail revolve: a disc whose rim is an ellipse with semi-axes 2 and 1
rail = nf.make_circle((0, 0, 0), (0, 0, 1), 1.0).transformed(
    nf.Transform.scaling((0, 0, 0), (2, 1, 1)))
disc = nf.rail_revolve(nf.make_line((0, 0, 0), (2, 0, 0)), rail, (0, 0, 0), (0, 0, 1))
print(f"rail revolve: area {nf.TrimmedSurface(disc).area():.6f} (ellipse pi a b = "
      f"{2 * math.pi:.6f})")

# --- tapered extrusion of a square: four planar faces, closed with caps -> frustum
square = nf.make_rectangle_curve((0, 0, 0), (2, 0, 0), (0, 2, 0))
angle = math.radians(10)
faces = nf.extrude_tapered(square, (0, 0, 1), angle)
frustum = nf.cap(nf.Shell.sew([nf.TrimmedSurface(f) for f in faces]))
a1, a2 = 4.0, (2 + 2 * math.tan(angle)) ** 2
print(f"tapered extrusion: {len(faces)} faces, volume {abs(frustum.volume()):.9f} "
      f"(h/3 (A1 + A2 + sqrt(A1 A2)) = {(a1 + a2 + math.sqrt(a1 * a2)) / 3:.9f})")

# --- ribbon along an L-shaped polyline: mitred at the corner
strip = nf.ribbon(nf.make_polyline([(0, 0, 0), (2, 0, 0), (2, 1, 0)]), 0.5, (0, 0, 1))
print(f"ribbon: {len(strip)} pieces, area {sum(nf.TrimmedSurface(s).area() for s in strip):.9f}"
      " (1.25)")

# --- cone over a circle (exact)
cone = nf.extrude_to_point(nf.make_circle((0, 0, 0), (0, 0, 1), 1.0), (0, 0, 2))
print(f"cone: area {nf.TrimmedSurface(cone).area():.9f} (pi r s = {math.pi * math.sqrt(5):.9f})")

# --- Coons patch from two arcs and two lines: an exact cylinder patch
a0 = nf.make_arc((0, 0, 0), (1, 0, 0), (0, 1, 0), 1.0, 0, math.pi / 2)
a1 = nf.make_arc((0, 0, 1), (1, 0, 0), (0, 1, 0), 1.0, 0, math.pi / 2)
patch = nf.coons(a0, a1, nf.make_line((1, 0, 0), (1, 0, 1)), nf.make_line((0, 1, 0), (0, 1, 1)))
dev = max(abs(math.hypot(p.x, p.y) - 1)
          for p in (patch(i / 10, j / 10) for i in range(11) for j in range(11)))
print(f"coons from arcs: rational={bool(patch.weights)}, max |r - 1| = {dev:.1e}")

# --- three-sided patch bounded by three quarter circles of the unit sphere
arcs = [nf.make_arc((0, 0, 0), x, y, 1.0, 0, math.pi / 2)
        for x, y in (((1, 0, 0), (0, 1, 0)), ((0, 1, 0), (0, 0, 1)), ((0, 0, 1), (1, 0, 0)))]
tri = nf.triangle_patch(*arcs)
edge = max(abs(math.sqrt(p.x ** 2 + p.y ** 2 + p.z ** 2) - 1)
           for p in (tri(t / 10, 0) for t in range(11)))
print(f"triangle patch: boundary on the sphere within {edge:.1e}, corner {tri(0.5, 1)}")

# --- rebuild: a circle as a uniform cubic B-spline
for n in (8, 16, 32):
    r = nf.rebuild_curve(nf.make_circle((0, 0, 0), (0, 0, 1), 1.0), n)
    print(f"rebuild circle with {n:2d} control points: max error {r.max_error:.2e}")
