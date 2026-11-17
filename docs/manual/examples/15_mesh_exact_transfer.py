# requires: gmsh
"""Exact NURBS transfer to gmsh: surfaces, trimmed surfaces, curves; 2D meshing."""
import math

import gmsh
import nurbforge as nf
from nurbforge.mesh import GmshSession

# a sphere of radius 2 as one rational surface (revolved semicircle)
meridian = nf.make_arc((0, 0, 0), (0, 0, -1), (1, 0, 0), 2.0, 0.0, math.pi)
sphere = nf.revolve(meridian, (0, 0, 0), (0, 0, 1), 2 * math.pi)

# a flat plate with a circular hole: trimmed surface, loops in (u, v)
plate = nf.TrimmedSurface(nf.make_rectangle((3, -1, 0), (2, 0, 0), (0, 2, 0)))
plate = plate.with_hole(nf.circle_loop(0.5, 0.5, 0.2))

with GmshSession() as g:
    s_tag = g.add_surface(sphere)        # NurbsSurface -> add_surface
    p_tag = g.add_surface(plate)         # TrimmedSurface -> add_trimmed_surface
    c_tag = g.add_curve(nf.make_circle((0, 0, 3), (0, 0, 1), 1.0))
    gmsh.model.occ.synchronize()

    # gmsh evaluates the transferred surface exactly like Nurbforge
    worst = 0.0
    for u in (0.1, 0.37, 0.8):
        for v in (0.2, 0.5, 0.9):
            p = sphere(u, v)
            q = gmsh.model.getValue(2, s_tag, [u, v])
            worst = max(worst, math.dist((p.x, p.y, p.z), q))
    print(f"max |gmsh - nurbforge| on the sphere: {worst:.2e}")

    g.mesh(dim=2, size=0.25)
    for tag, name in [(s_tag, "sphere"), (p_tag, "plate")]:
        types, elems, _ = gmsh.model.mesh.getElements(2, tag)
        print(f"{name}: {len(elems[0])} triangles")
    nodes, xyz, _ = gmsh.model.mesh.getNodes(2, s_tag, includeBoundary=True)
    r = [math.hypot(*xyz[3 * i:3 * i + 3]) for i in range(len(nodes))]
    print(f"sphere nodes: radius between {min(r):.12f} and {max(r):.12f}")
    g.write("transfer.msh")
    g.write("transfer.vtk")
print("wrote transfer.msh and transfer.vtk")
