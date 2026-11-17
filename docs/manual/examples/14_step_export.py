"""Write B-rep shells to STEP with the native writer (no gmsh, no OpenCASCADE)."""
import re

import nurbforge as nf

# a box with a bore: a closed shell with trimmed planar faces and a cylindrical face
block = nf.make_box((-2, -2, 0), (4, 4, 2))
bore = nf.make_cylinder((0, 0, -1), (0, 0, 1), 0.5, 4)
part = nf.boolean(block, bore, nf.BooleanOp.Difference)
print(f"part: {len(part.faces)} faces, closed={part.is_closed}, volume={part.volume():.6f}")

# write a file; the timestamp is fixed so that the file is reproducible
options = {"unit": "mm", "name": "block", "author": "manual", "organization": "Nurbforge",
           "timestamp": "2026-01-01T00:00:00", "uncertainty": 1e-6}
nf.write_step("block.step", part, **options)

# the same contents as a string, e.g. to store in a database or to inspect
text = nf.to_step(part, **options)
with open("block.step", encoding="utf-8") as f:
    print("file and string identical:", f.read() == text)

# count the main entities of the file
for entity in ["MANIFOLD_SOLID_BREP", "CLOSED_SHELL", "ADVANCED_FACE", "EDGE_CURVE",
               "VERTEX_POINT", "B_SPLINE_SURFACE_WITH_KNOTS", "B_SPLINE_CURVE_WITH_KNOTS"]:
    n = len(re.findall(r"\b" + entity + r"\(", text))
    print(f"{entity:28s} {n}")
print("schema:", re.search(r"FILE_SCHEMA\(\('([^']*)'", text).group(1))

# several shells, each with its own solid name; an open shell becomes a surface model
lid = nf.make_cylinder((0, 0, 2), (0, 0, 1), 0.5, 0.2)
open_shell = nf.Shell.sew(list(nf.make_box((5, 0, 0), (1, 1, 1)).faces)[:-1])
nf.write_step("assembly.step", [part, lid, open_shell], unit="mm",
              names=["Block", "Plug", "Cover"])
print("open shell closed:", open_shell.is_closed)
for name, shell in nf.read_step_named("assembly.step"):
    print(f"  {name:6s} faces={len(shell.faces):2d} closed={shell.is_closed}")
