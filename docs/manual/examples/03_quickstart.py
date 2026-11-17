"""Quick start: from a curve to a drilled solid, its mass properties and a STEP file."""
import math

import nurbforge as nf

# 1. A curve: a quarter circle of radius 10 in the xz plane (exact rational NURBS)
arc = nf.make_arc(center=(0, 0, 0), x_axis=(1, 0, 0), y_axis=(0, 0, 1),
                  radius=10.0, start_angle=0.0, end_angle=math.pi / 2)
print(f"arc: degree {arc.degree}, {len(arc.control_points)} control points, "
      f"rational {arc.is_rational}, domain {arc.domain}")
print(f"arc length {nf.length(arc):.12f} (exact {5 * math.pi:.12f})")
mid = arc(0.5)
print(f"arc(0.5) = ({mid.x:.6f}, {mid.y:.6f}, {mid.z:.6f}), |arc(0.5)| = {nf.norm(mid):.15f}")

# 2. A surface: revolving the arc about z gives a hemisphere
dome = nf.revolve(arc, axis_point=(0, 0, 0), axis_dir=(0, 0, 1), angle=2 * math.pi)
print(f"dome: degrees ({dome.degree_u}, {dome.degree_v}), net {dome.num_u} x {dome.num_v}")
print(f"dome area {nf.area(dome):.10f} (exact {200 * math.pi:.10f})")

# 3. Solids: a box with a cylindrical bore (Boolean difference)
box = nf.make_box(corner=(-20, -20, 0), size=(40, 40, 10))
bore = nf.make_cylinder(base=(0, 0, -5), axis=(0, 0, 1), radius=6.0, height=20.0)
part = nf.boolean(box, bore, nf.BooleanOp.Difference)
print(f"part: {len(part.faces)} faces, closed {part.is_closed}, labels {sorted(part.face_labels)}")

# 4. Mass properties (unit density)
mp = part.mass_properties()
exact = 40 * 40 * 10 - math.pi * 6.0**2 * 10
print(f"volume {mp.volume:.8f} (exact {exact:.8f})")
c = mp.centroid
print(f"centroid ({c.x:.6f}, {c.y:.6f}, {c.z:.6f})")
print(f"inside at (0, 0, 5): {nf.point_in_shell(part, (0, 0, 5))}, "
      f"at (10, 10, 5): {nf.point_in_shell(part, (10, 10, 5))}")

# 5. STEP export with the native writer, and reading it back
nf.write_step("part.step", part, unit="mm", name="drilled block")
back = nf.read_step("part.step")
print(f"read back: {len(back)} shell(s), volume {back[0].volume():.8f}")
