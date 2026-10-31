"""Biconvex lens surfaces built by revolution, meshed with gmsh.

Requires: pip install gmsh
"""
import math

import nurbforge as nf
from nurbforge.mesh import mesh_surfaces

R1, R2 = 50.0, 80.0          # radii of curvature [mm]
aperture = 25.0              # clear diameter [mm]
center_thickness = 5.0
h = aperture / 2

def spherical_cap(R, z_vertex, sign):
    """Cap with vertex at z_vertex; sign=+1 curves towards -z, -1 towards +z."""
    alpha = math.asin(h / R)
    center = (0, 0, z_vertex - sign * R)
    meridian = nf.make_arc(center, (0, 0, sign), (1, 0, 0), R, 0.0, alpha)
    return nf.revolve(meridian, (0, 0, 0), (0, 0, 1), 2 * math.pi)

front = spherical_cap(R1, center_thickness / 2, +1)
back = spherical_cap(R2, -center_thickness / 2, -1)
edge_z0 = back(0.0, 1.0).z
edge_z1 = front(0.0, 1.0).z
edge = nf.revolve(nf.make_line((h, 0, edge_z0), (h, 0, edge_z1)), (0, 0, 0), (0, 0, 1), 2 * math.pi)

mesh_surfaces([front, back, edge], "lens.msh", size=1.0)
print("wrote lens.msh")
