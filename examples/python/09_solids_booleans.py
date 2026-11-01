"""Solids: sew faces into closed shells, volumes, point-in-solid and Boolean operations.

Optional: writes drilled_cube.step if gmsh is installed.
"""
import math

import nurbforge as nf


def box(h):
    """Cube [-h, h]^3 from six planar faces (Shell.oriented() fixes the face orientations)."""
    s = 2 * h
    rects = [((-h, -h, -h), (s, 0, 0), (0, s, 0)), ((-h, -h, h), (s, 0, 0), (0, s, 0)),
             ((-h, -h, -h), (s, 0, 0), (0, 0, s)), ((-h, h, -h), (s, 0, 0), (0, 0, s)),
             ((-h, -h, -h), (0, s, 0), (0, 0, s)), ((h, -h, -h), (0, s, 0), (0, 0, s))]
    return nf.Shell.sew([nf.TrimmedSurface(nf.make_rectangle(*r)) for r in rects]).oriented()


def ball(center, r):
    """Sphere from a revolved semicircle: one face, closed by its seam edge and two poles."""
    meridian = nf.make_arc(center, (0, 0, -1), (1, 0, 0), r, 0.0, math.pi)
    return nf.Shell.sew([nf.TrimmedSurface(nf.revolve(meridian, center, (0, 0, 1), 2 * math.pi))])


def cylinder(r, z0, z1):
    """Cylinder along z: extruded circle plus two lids trimmed to circles."""
    side = nf.extrude(nf.make_circle((0, 0, z0), (0, 0, 1), r), (0, 0, z1 - z0))
    lids = [nf.TrimmedSurface(nf.make_rectangle((-r, -r, z), (2 * r, 0, 0), (0, 2 * r, 0)),
                              nf.circle_loop(0.5, 0.5, 0.5)) for z in (z0, z1)]
    return nf.Shell.sew([nf.TrimmedSurface(side), *lids]).oriented()


cube = box(1.0)
print(f"cube: {len(cube.faces)} faces, {len(cube.edges)} edges, {len(cube.vertices)} vertices, "
      f"closed={cube.is_closed}, oriented={cube.is_consistently_oriented}, "
      f"volume {cube.volume():.12f}")
print(f"point_in_shell(cube, (0.9, 0, 0)) = {nf.point_in_shell(cube, (0.9, 0, 0))}, "
      f"(1.1, 0, 0) = {nf.point_in_shell(cube, (1.1, 0, 0))}")

# --- two unit spheres at distance d: union, intersection (lens), difference ----------------
R, d = 1.0, 1.0
a, b = ball((0, 0, 0), R), ball((d, 0, 0), R)
sphere_vol = 4 / 3 * math.pi * R**3
lens = math.pi * (4 * R + d) * (2 * R - d) ** 2 / 12  # analytic lens volume
expected = {nf.BooleanOp.Union: 2 * sphere_vol - lens, nf.BooleanOp.Intersection: lens,
            nf.BooleanOp.Difference: sphere_vol - lens}
for op, ref in expected.items():
    s = nf.boolean(a, b, op)
    vol = s.volume(1e-7)  # relative tolerance of the volume integral
    print(f"spheres {op.name:<12}: {len(s.faces)} faces, closed={s.is_closed}, "
          f"volume {vol:.8f} (analytic {ref:.8f}, rel. error {abs(vol - ref) / ref:.1e})")

# --- drill a cylinder r = 0.3 through the cube ----------------------------------------------
r = 0.3
drilled = nf.boolean(cube, cylinder(r, -2.0, 2.0), nf.BooleanOp.Difference)
ref = 8.0 - math.pi * r * r * 2.0
print(f"drilled cube: {len(drilled.faces)} faces, volume {drilled.volume(1e-7):.8f} "
      f"(8 - 2 pi r^2 = {ref:.8f}), origin inside: "
      f"{nf.point_in_shell(drilled, (0, 0, 0))}")

try:  # optional STEP export through gmsh's OpenCASCADE kernel
    from nurbforge.mesh import export_step

    export_step([drilled], "drilled_cube.step")
    print("wrote drilled_cube.step")
except (ImportError, OSError) as exc:  # gmsh missing or its runtime libraries unavailable
    print(f"STEP export skipped: {exc}")
