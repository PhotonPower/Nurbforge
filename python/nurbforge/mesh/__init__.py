"""Meshing via gmsh (optional dependency: ``pip install nurbforge[mesh]``)."""
from .gmsh_export import (  # noqa: F401
    GmshSession,
    add_curve,
    add_shell,
    add_surface,
    add_trimmed_surface,
    export_step,
    mesh_surfaces,
)
