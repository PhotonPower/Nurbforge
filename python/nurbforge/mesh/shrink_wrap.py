"""Shrink wrap: one closed triangle mesh around arbitrary input (approximate, mesh level).

:func:`shrink_wrap` builds a watertight envelope around shells and loose surfaces, whether they
overlap, leave gaps or have naked edges. Unlike everything else in Nurbforge the result is
**not exact geometry**: it is a triangle mesh accurate to about the ``resolution``.

The work is done by the C++ kernel (``nf::shrink_wrap``, no dependencies): the input is sampled
into points, the unsigned distance to them is computed on a regular grid in a narrow band, the
grid nodes that can be reached from outside without coming closer to the input than the
envelope distance are marked as outside, and the surface between outside and the rest is
extracted with marching tetrahedra, which always gives a closed, consistently oriented
2-manifold.
"""
from __future__ import annotations

from collections.abc import Iterable

import numpy as np

from .._core import NurbsSurface, Shell, TrimmedSurface, _shrink_wrap
from .surface_mesh import SurfaceMesh


def _faces(objects) -> list[TrimmedSurface]:
    if isinstance(objects, (Shell, TrimmedSurface, NurbsSurface)):
        objects = [objects]
    elif not isinstance(objects, Iterable):
        raise TypeError("shrink_wrap: expected shells or surfaces")
    faces: list[TrimmedSurface] = []
    for obj in objects:
        if isinstance(obj, Shell):
            faces.extend(obj.faces)
        elif isinstance(obj, TrimmedSurface):
            faces.append(obj)
        elif isinstance(obj, NurbsSurface):
            faces.append(TrimmedSurface(obj))
        else:
            raise TypeError(f"shrink_wrap: cannot wrap a {type(obj).__name__}; "
                            "expected Shell, TrimmedSurface or NurbsSurface")
    return faces


def shrink_wrap(objects, resolution: float, offset: float = 0.0, fill_holes: bool = True,
                max_nodes: int = 20_000_000, min_quality: float = 0.3) -> SurfaceMesh:
    """One closed triangle mesh around shells and surfaces (approximate).

    ``objects``: a ``Shell``, ``TrimmedSurface`` or ``NurbsSurface``, or a list of them. They
    may overlap, leave gaps or be open. ``resolution`` is the cell size of the grid the envelope
    is built on, which is about the edge length of the mesh and the accuracy of the result:
    openings and gaps narrower than about 1.5 x ``resolution`` (plus twice the ``offset``) are
    closed, smaller details are smoothed away. ``offset`` is the distance of the envelope from
    the input.

    The envelope is first built at a distance of at least 0.75 x ``resolution``, which is what
    makes it watertight on the grid; for a smaller ``offset`` the surface that lies the
    difference further in is taken (a morphological closing). Parts thinner than about 1.5 x
    ``resolution`` do not survive this: an open surface needs an ``offset`` of at least
    0.75 x ``resolution``, a ``ValueError`` is raised if nothing is left.

    With ``fill_holes`` (default) everything that cannot be reached from outside is filled,
    and the mesh is the outer envelope only; without it, enclosed voids larger than the
    envelope distance are kept as inner components (normals pointing into the void). Note that
    the inside of a closed shell is such a void: the input is a set of surfaces.

    The triangles cut out of the grid can be arbitrarily thin. If the worst of them has a
    :meth:`SurfaceMesh.quality` below ``min_quality``, the nodes are moved within the surface
    until it is reached (at most 50 smoothing steps; typically 5 give a worst quality above
    0.3 and a mean of 0.85). ``min_quality=0`` returns the mesh as cut out of the grid. For
    triangles of a different size see :func:`nurbforge.mesh.remesh`.

    Returns a :class:`SurfaceMesh` with linear triangles, watertight and with normals pointing
    out of the wrapped volume. Parts further apart than the gap size give separate closed
    components. Raises ``ValueError`` if the grid would exceed ``max_nodes`` nodes.
    """
    nodes, triangles = _shrink_wrap(_faces(objects), float(resolution), float(offset),
                                    bool(fill_holes), int(max_nodes), float(min_quality))
    return SurfaceMesh(nodes, triangles, np.zeros(len(triangles), dtype=np.int64),
                       np.arange(1, len(nodes) + 1, dtype=np.int64))
