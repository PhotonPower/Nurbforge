"""Meshing via gmsh (optional dependency: ``pip install nurbforge[mesh]``) and mesh-level tools
that need numpy only (``shrink_wrap``, ``SurfaceMesh``)."""
from .gmsh_export import (  # noqa: F401
    GmshSession,
    add_curve,
    add_shell,
    add_surface,
    add_trimmed_surface,
    bem_mesh,
    export_step,
    mesh_surfaces,
    remesh,
)
from .shrink_wrap import shrink_wrap  # noqa: F401
from .surface_mesh import MeshTopologyError, SurfaceMesh  # noqa: F401
