# requires: gmsh
"""A biconvex asphere lens: two optical surfaces and the cylindrical edge, meshed with gmsh."""

import nurbforge as nf
from nurbforge.mesh import GmshSession

r_max, centre_thickness = 12.0, 6.0
front_sag = nf.Asphere(curvature=1 / 40.0, conic=-0.8, even=[1e-6, -2e-9])
back_sag = nf.Asphere(curvature=-1 / 120.0)  # R < 0: centre on -z, the back bulges to +z

front = nf.asphere_surface(front_sag, r_max, tolerance=1e-8)  # vertex at z = 0, opens to +z
back = nf.asphere_surface(back_sag, r_max).transformed(
    nf.Transform.translation((0, 0, centre_thickness)))
# the edge: a cylinder between the two rims
z0 = front_sag.sag(r_max)
z1 = centre_thickness + back_sag.sag(r_max)
edge = nf.extrude(nf.make_circle((0, 0, z0), (0, 0, 1), r_max), (0, 0, z1 - z0))
print(f"rim heights: front {z0:.6f}, back {z1:.6f}, edge thickness {z1 - z0:.6f}")
print(f"areas: front {nf.area(front):.6f}, back {nf.area(back):.6f}, edge {nf.area(edge):.6f}")

with GmshSession() as g:
    for s in (front, back, edge):
        g.add_surface(s)
    g.mesh(dim=2, size=1.0)
    g.write("lens.msh")
    _, coords, _ = g.gmsh.model.mesh.getNodes()
print(f"wrote lens.msh with {len(coords) // 3} nodes")
