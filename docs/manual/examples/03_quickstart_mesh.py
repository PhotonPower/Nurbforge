# requires: gmsh
"""Quick start, continued: mesh the drilled block with gmsh (volume and BEM surface mesh)."""
import nurbforge as nf
from nurbforge.mesh import GmshSession, bem_mesh

box = nf.make_box((-20, -20, 0), (40, 40, 10))
bore = nf.make_cylinder((0, 0, -5), (0, 0, 1), 6.0, 20.0)
part = nf.boolean(box, bore, nf.BooleanOp.Difference)

# tetrahedral volume mesh with one physical group per face label
with GmshSession() as g:
    g.add_shell(part, physical_groups=True)
    g.refine("side", 1.0)            # finer elements on the bore wall
    g.mesh(dim=3, size=4.0)
    g.write("part.msh")
    print("face labels -> gmsh surfaces:", g.face_tags())

# closed, outward-oriented triangle mesh for a boundary element solver
m = bem_mesh(part, size=3.0)
print(f"BEM mesh: {len(m.nodes)} nodes, {len(m.triangles)} triangles, "
      f"watertight {m.is_watertight()}, volume {m.volume():.3f} (exact {part.volume():.3f})")
