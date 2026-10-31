"""Meshing via gmsh (optional dependency: ``pip install nurbforge[mesh]``)."""
from .gmsh_export import GmshSession, add_curve, add_surface, mesh_surfaces  # noqa: F401
