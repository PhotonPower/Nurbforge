# requires: gmsh
"""BEM-ready surface meshes: watertight, outward-oriented, as numpy arrays."""
import math

import numpy as np
import nurbforge as nf
from nurbforge.mesh import GmshSession, MeshTopologyError, bem_mesh

# a box with a spherical cavity: two closed components
part = nf.boolean(nf.make_box((-2, -2, -2), (4, 4, 4)), nf.make_sphere((0, 0, 0), 1.0),
                  nf.BooleanOp.Difference)
m = bem_mesh(part, size=0.3, order=2, refine={"surface": 0.15}, path="cavity.msh")
print(f"nodes {m.nodes.shape}, triangles {m.triangles.shape}, order {m.order}")
print("groups:", {k: len(v) for k, v in m.groups.items()})
print(f"watertight {m.is_watertight()}, consistent {m.is_consistently_oriented()}")
print(f"mesh volume {m.volume():.4f}, exact {64 - 4 / 3 * math.pi:.4f}")

# normals point out of the solid: away from the box, into the cavity
c = m.centroids()
radial = np.einsum("ij,ij->i", m.normals(), c)
inner = np.linalg.norm(c, axis=1) < 1.5
print(f"outer faces pointing outwards: {np.all(radial[~inner] > 0)}")
print(f"cavity faces pointing to the centre: {np.all(radial[inner] < 0)}")
print(f"total area {m.areas().sum():.4f}, exact {6 * 16 + 4 * math.pi:.4f}")

# the same step by step in a session; orient=False shows what gmsh produced
with GmshSession() as g:
    g.add_shell(nf.make_cylinder((0, 0, 0), (0, 0, 1), 1.0, 2.0), solid=False)
    g.mesh(dim=2, size=0.2)
    g.gmsh.model.mesh.reverse([(2, g.face_tags()["top"][0])])  # spoil one face
    print("before:", g.surface_mesh(orient=False).is_consistently_oriented())
    mesh = g.surface_mesh()                                     # repairs, also in gmsh
    print("after:", mesh.is_consistently_oriented(), f"volume {mesh.volume():.4f}")
    g.write("cylinder.stl")

# an open shell cannot be a BEM boundary
box = nf.make_box((0, 0, 0), (1, 1, 1))
try:
    bem_mesh(nf.Shell.sew(list(box.faces)[:-1]), size=0.3)
except MeshTopologyError as exc:
    print("MeshTopologyError:", exc)
