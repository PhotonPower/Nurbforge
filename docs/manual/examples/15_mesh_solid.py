# requires: gmsh
"""Shells, physical groups from face labels, local mesh sizes and a volume mesh."""
import gmsh
import numpy as np
import nurbforge as nf
from nurbforge.mesh import GmshSession

# a pipe with labelled faces and a slot cut by a Boolean: the labels survive the Boolean
pipe = nf.make_cylinder((0, 0, 0), (0, 0, 1), 1.0, 5.0).with_face_labels(
    ["wall", "inlet", "outlet"])
slot = nf.make_box((-2, -0.2, 2), (4, 0.4, 1))
part = nf.boolean(pipe, slot, nf.BooleanOp.Difference)
print("face labels:", sorted(set(part.face_labels)))

with GmshSession() as g:
    tags = g.add_shell(part, physical_groups=True, volume_name="fluid")
    print(f"surfaces {tags[:-1]}, volume {tags[-1]}")
    print("face_tags:", g.face_tags())
    g.refine("inlet", 0.05, distance=1.0)   # fine at the inlet, growing over 1.0
    g.refine("outlet", 0.1)                  # only the outlet face itself
    g.mesh(dim=3, size=0.3, curvature=24)    # tetrahedra, 24 elements per 2 pi
    for dim, tag in gmsh.model.getPhysicalGroups():
        name = gmsh.model.getPhysicalName(dim, tag)
        print(f"  physical group dim={dim} {name!r}")

    # volume of the tetrahedral mesh against the exact volume
    _, _, conn = gmsh.model.mesh.getElements(3)
    ids, xyz, _ = gmsh.model.mesh.getNodes()
    pos = dict(zip(ids, np.asarray(xyz).reshape(-1, 3)))
    tets = np.array([pos[n] for n in conn[0]]).reshape(-1, 4, 3)
    a, b, c, d = (tets[:, i] for i in range(4))
    vol = np.abs(np.einsum("ij,ij->i", b - a, np.cross(c - a, d - a))).sum() / 6
    print(f"{len(tets)} tetrahedra, mesh volume {vol:.4f}, exact {part.volume():.4f}")
    g.write("pipe.msh")
print("wrote pipe.msh")
