"""Solid primitives, Boolean operations, point classification, face labels and voids."""
import math

import nurbforge as nf

# --- primitives: closed shells with exact faces, outward normals and labelled
prims = {
    "box": nf.make_box((0, 0, 0), (2, 3, 4)),
    "sphere": nf.make_sphere((0, 0, 0), 1.0),
    "cylinder": nf.make_cylinder((0, 0, 0), (0, 0, 1), 1.0, 2.0),
    "cone": nf.make_cone((0, 0, 0), (0, 0, 1), 1.0, 0.0, 3.0),
    "frustum": nf.make_cone((0, 0, 0), (0, 0, 1), 1.0, 0.5, 2.0),
    "torus": nf.make_torus((0, 0, 0), (0, 0, 1), 2.0, 0.5),
}
exact = {"box": 24, "sphere": 4 / 3 * math.pi, "cylinder": 2 * math.pi, "cone": math.pi,
         "frustum": math.pi * 2 / 3 * (1 + 0.5 + 0.25), "torus": 2 * math.pi**2 * 2 * 0.25}
for name, s in prims.items():
    print(f"{name:8s}: {len(s.faces)} faces {s.face_labels}, volume {s.volume():.10f} "
          f"(exact {exact[name]:.10f})")

# --- Booleans ---------------------------------------------------------------
block = nf.make_box((-2, -2, 0), (4, 4, 2))
bore = nf.make_cylinder((0, 0, -1), (0, 0, 1), 0.5, 4)
drilled = nf.boolean(block, bore, nf.BooleanOp.Difference)
print(f"block - bore: {len(drilled.faces)} faces, closed={drilled.is_closed}, "
      f"volume {drilled.volume(1e-9):.10f} (32 - 2 pi 0.25 = {32 - 0.5 * math.pi:.10f})")
print(f"labels after the Boolean: {drilled.face_labels}")
print(f"faces labelled 'side' (the bore wall): {drilled.faces_labeled('side')}")

a = nf.make_sphere((0, 0, 0), 1.0)
b = nf.make_sphere((1, 0, 0), 1.0)
lens = math.pi * (4 + 1) * (2 - 1) ** 2 / 12
for op, ref in ((nf.BooleanOp.Union, 8 / 3 * math.pi - lens), (nf.BooleanOp.Intersection, lens),
                (nf.BooleanOp.Difference, 4 / 3 * math.pi - lens)):
    s = nf.boolean(a, b, op)
    v = s.volume(1e-9)
    print(f"spheres {op.name:12s}: volume {v:.9f}, exact {ref:.9f}, rel. error {abs(v - ref) / ref:.1e}")

# --- point classification ---------------------------------------------------
for p in ((0, 0, 1), (1.5, 1.5, 1), (0, 0, 3)):
    print(f"point_in_shell(drilled, {p}) = {nf.point_in_shell(drilled, p)}")

# --- labels name faces for boundary conditions; they survive the operations -
pipe = nf.make_cylinder((0, 0, 0), (0, 0, 1), 1.0, 5.0).with_face_labels(["wall", "inlet", "outlet"])
slot = nf.make_box((-2, -0.2, 2), (4, 0.4, 1)).with_face_labels(["slot"] * 6)
part = nf.boolean(pipe, slot, nf.BooleanOp.Difference)
part = part.transformed(nf.Transform.translation((10, 0, 0)))
groups = {lab: len(part.faces_labeled(lab)) for lab in dict.fromkeys(part.face_labels)}
print(f"slotted pipe: {len(part.faces)} faces, faces per label {groups}")

# --- a hollow body: box minus an inner sphere -> two components, the void poi
hollow = nf.boolean(nf.make_box((-2, -2, -2), (4, 4, 4)), nf.make_sphere((0, 0, 0), 1.0),
                    nf.BooleanOp.Difference)
comps = hollow.components()
print(f"hollow box: {len(comps)} components, volumes "
      f"{', '.join(f'{c.volume():+.6f}' for c in comps)}, total {hollow.volume():.6f} "
      f"(64 - 4/3 pi = {64 - 4 / 3 * math.pi:.6f})")
print(f"centre inside the solid? {nf.point_in_shell(hollow, (0, 0, 0))}")
