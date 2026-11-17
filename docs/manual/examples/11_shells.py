"""B-rep shells: sewing faces, topology queries, orientation, volume and mass properties."""
import math

import nurbforge as nf

# --- a cylinder sewn from three faces ---------------------------------------
r, h = 1.0, 2.0
side = nf.TrimmedSurface(nf.extrude(nf.make_circle((0, 0, 0), (0, 0, 1), r), (0, 0, h)))
lids = [nf.TrimmedSurface(nf.make_rectangle((-r, -r, z), (2 * r, 0, 0), (0, 2 * r, 0)),
                          nf.circle_loop(0.5, 0.5, 0.5)) for z in (0.0, h)]
raw = nf.Shell.sew([side, *lids], tolerance=1e-6)
print(f"sewn: {len(raw.faces)} faces, {len(raw.edges)} edges, {len(raw.vertices)} vertices, "
      f"closed={raw.is_closed}, consistently oriented={raw.is_consistently_oriented}")
cyl = raw.oriented()  # flips faces until all normals point outwards
print(f"oriented: consistently oriented={cyl.is_consistently_oriented}, "
      f"volume {cyl.volume():.10f} (pi r^2 h = {math.pi * r * r * h:.10f})")

# --- topology: edges and their uses -----------------------------------------
for i, e in enumerate(cyl.edges):
    uses = ", ".join(f"face {u.face} loop {u.loop}{' rev' if u.reversed else ''}" for u in e.uses)
    kind = "closed" if e.v0 == e.v1 else "open"
    print(f"edge {i}: {kind}, length {nf.length(e.curve):.6f}, used by [{uses}]")

# --- an open shell: leave out the top lid -----------------------------------
cup = nf.Shell.sew([side, lids[0]])
print(f"cup: closed={cup.is_closed}, free edges {cup.free_edges()}, "
      f"non-manifold edges {cup.non_manifold_edges()}, area {cup.area():.8f}")
try:
    cup.volume()
except nf.GeometryError as exc:
    print(f"volume of an open shell: GeometryError: {exc}")

# --- mass properties --------------------------------------------------------
mp = cyl.mass_properties()
print(f"volume {mp.volume:.8f}, centroid ({mp.centroid.x:.6f}, {mp.centroid.y:.6f}, "
      f"{mp.centroid.z:.6f})")
print("inertia about the centroid (unit density):")
for row in mp.inertia:
    print("   " + "  ".join(f"{x:10.6f}" for x in row))
m = math.pi * r * r * h
print(f"analytic: Ixx = m (3 r^2 + h^2) / 12 = {m * (3 * r * r + h * h) / 12:.6f}, "
      f"Izz = m r^2 / 2 = {m * r * r / 2:.6f}")
print(f"principal moments {tuple(round(x, 6) for x in mp.principal_moments)}")

# --- transforms keep the orientation (mirrors flip faces back) --------------
moved = cyl.transformed(nf.Transform.rotation((0, 0, 0), (1, 0, 0), math.pi / 2))
mirrored = cyl.transformed(nf.Transform.mirror((0, 0, 0), (1, 0, 0)))
print(f"rotated volume {moved.volume():.8f}, mirrored volume {mirrored.volume():.8f}")
