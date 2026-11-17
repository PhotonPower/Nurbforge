"""Read STEP files with the native reader: units, tolerance, names."""
import warnings

import nurbforge as nf

# a 2 x 2 x 2 cube with a spherical void, written in metres
cube = nf.make_box((-1, -1, -1), (2, 2, 2))
void = nf.make_sphere((0, 0, 0), 0.5)
hollow = nf.boolean(cube, void, nf.BooleanOp.Difference)
print(f"hollow: components={len(hollow.components())}, volume={hollow.volume():.9f}")
nf.write_step("hollow_m.step", hollow, unit="m")

# read back in different target units: coordinates are converted
for unit, scale in [("m", 1.0), ("mm", 1e3), ("cm", 1e2)]:
    (shell,) = nf.read_step("hollow_m.step", unit=unit)
    print(f"read as {unit:2s}: volume = {shell.volume():.6g} "
          f"(expected {hollow.volume() * scale**3:.6g}), closed={shell.is_closed}")

# names survive a round trip in micrometres (e.g. parts for two-photon lithography)
pillar = nf.make_cylinder((0, 0, 0), (0, 0, 1), 5.0, 20.0)
cap = nf.make_sphere((0, 0, 25), 6.0)
text = nf.to_step([pillar, cap], unit="um", names=["Pillar", "Cap"])
for name, shell in nf.from_step_named(text, unit="um"):
    print(f"{name:7s} volume = {shell.volume():.4f} um^3")

# an unsupported entity raises GeometryError naming entity, instance and face ...
box_text = nf.to_step(nf.make_box((0, 0, 0), (1, 1, 1)))
i = box_text.index("B_SPLINE_SURFACE_WITH_KNOTS(")
broken = box_text[:i] + "OFFSET_SURFACE(" + box_text[i + len("B_SPLINE_SURFACE_WITH_KNOTS("):]
try:
    nf.from_step(broken)
except nf.GeometryError as exc:
    print("GeometryError:", exc)

# ... or, with skip_unsupported, the face is left out with a UserWarning
with warnings.catch_warnings(record=True) as caught:
    warnings.simplefilter("always")
    (shell,) = nf.from_step(broken, skip_unsupported=True)
print(f"skipped: {len(caught)} warning(s), {len(shell.faces)} faces, closed={shell.is_closed}")
