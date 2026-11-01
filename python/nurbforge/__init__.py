"""Nurbforge - headless NURBS-based CAD kernel for scientific software.

The geometry kernel is implemented in C++ (``nurbforge._core``); this package adds
Pythonic helpers, meshing via gmsh (:mod:`nurbforge.mesh`) and an optional LLM
interface via Ollama (:mod:`nurbforge.llm`).
"""
from __future__ import annotations

from ._core import (  # noqa: F401
    CurveEnd,
    CurvePoint,
    Extension,
    GeometryError,
    NurbsCurve,
    NurbsSurface,
    Parametrization,
    SurfaceCurvature,
    SurfacePoint,
    SurfaceSide,
    Tolerance,
    Transform,
    Vec3,
    __version__,
    approximate_curve,
    approximate_curve_tol,
    approximate_surface,
    arc_length,
    area,
    clamped_uniform_knots,
    closest_point,
    closest_point_local,
    cross,
    curvature,
    data_parameters,
    distance,
    dot,
    extend,
    extrude,
    grid_pattern,
    interpolate_curve,
    interpolate_surface,
    join,
    join_u,
    join_v,
    length,
    linear_pattern,
    loft,
    loft_parameters,
    make_arc,
    make_bilinear,
    make_circle,
    make_compatible,
    make_line,
    make_polyline,
    make_rectangle,
    norm,
    parameter_at_length,
    polar_pattern,
    reparametrize,
    revolve,
    ruled,
    surface_curvature,
    torsion,
    unique_knots,
    validate_knots,
)
from .patterns import apply_pattern  # noqa: F401

__all__ = [name for name in dir() if not name.startswith("_")]
