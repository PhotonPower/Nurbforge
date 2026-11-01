"""Whitelisted geometry tools exposed to the LLM (JSON-schema function calling)."""
from __future__ import annotations

import math
from collections.abc import Callable
from dataclasses import dataclass, field
from typing import Any

import numpy as np

from .. import _core as nf

Vec = list[float]


@dataclass
class Workspace:
    """Named object store the LLM operates on."""

    objects: dict[str, Any] = field(default_factory=dict)

    def put(self, name: str, obj: Any) -> str:
        self.objects[name] = obj
        return f"stored {name}: {obj!r}"

    def get(self, name: str) -> Any:
        if name not in self.objects:
            raise KeyError(f"unknown object '{name}'. Known: {sorted(self.objects)}")
        return self.objects[name]


@dataclass
class Tool:
    name: str
    description: str
    parameters: dict[str, Any]
    func: Callable[..., str]

    def schema(self) -> dict[str, Any]:
        return {"type": "function", "function": {
            "name": self.name, "description": self.description, "parameters": self.parameters}}


_VEC = {"type": "array", "items": {"type": "number"}, "minItems": 3, "maxItems": 3}
_NAME = {"type": "string", "description": "name of the object in the workspace"}


def _obj(props: dict[str, Any], required: list[str]) -> dict[str, Any]:
    return {"type": "object", "properties": props, "required": required}


def _deg(a: float) -> float:
    return math.radians(a)


def _build_tools() -> list[Tool]:
    T: list[Tool] = []

    def add(name, desc, props, req):
        def deco(f):
            T.append(Tool(name, desc, _obj(props, req), f))
            return f
        return deco

    @add("make_line", "Create a straight line segment.",
         {"name": _NAME, "p0": _VEC, "p1": _VEC}, ["name", "p0", "p1"])
    def _(ws, name, p0, p1):
        return ws.put(name, nf.make_line(p0, p1))

    @add("make_circle", "Create a full circle (exact NURBS).",
         {"name": _NAME, "center": _VEC, "normal": _VEC, "radius": {"type": "number"}},
         ["name", "center", "normal", "radius"])
    def _(ws, name, center, normal, radius):
        return ws.put(name, nf.make_circle(center, normal, radius))

    @add("make_arc", "Create a circular arc. Angles in degrees, measured from x_axis to y_axis.",
         {"name": _NAME, "center": _VEC, "x_axis": _VEC, "y_axis": _VEC,
          "radius": {"type": "number"}, "start_deg": {"type": "number"}, "end_deg": {"type": "number"}},
         ["name", "center", "x_axis", "y_axis", "radius", "start_deg", "end_deg"])
    def _(ws, name, center, x_axis, y_axis, radius, start_deg, end_deg):
        return ws.put(name, nf.make_arc(center, x_axis, y_axis, radius, _deg(start_deg), _deg(end_deg)))

    @add("make_rectangle", "Create a planar rectangular surface.",
         {"name": _NAME, "origin": _VEC, "edge_u": _VEC, "edge_v": _VEC},
         ["name", "origin", "edge_u", "edge_v"])
    def _(ws, name, origin, edge_u, edge_v):
        return ws.put(name, nf.make_rectangle(origin, edge_u, edge_v))

    @add("revolve", "Revolve a curve about an axis to create a surface. Angle in degrees.",
         {"name": _NAME, "profile": _NAME, "axis_point": _VEC, "axis_dir": _VEC,
          "angle_deg": {"type": "number"}},
         ["name", "profile", "axis_point", "axis_dir", "angle_deg"])
    def _(ws, name, profile, axis_point, axis_dir, angle_deg):
        return ws.put(name, nf.revolve(ws.get(profile), axis_point, axis_dir, _deg(angle_deg)))

    @add("extrude", "Extrude a curve along a direction vector to create a surface.",
         {"name": _NAME, "profile": _NAME, "direction": _VEC}, ["name", "profile", "direction"])
    def _(ws, name, profile, direction):
        return ws.put(name, nf.extrude(ws.get(profile), direction))

    @add("ruled", "Create a ruled surface between two curves.",
         {"name": _NAME, "curve0": _NAME, "curve1": _NAME}, ["name", "curve0", "curve1"])
    def _(ws, name, curve0, curve1):
        return ws.put(name, nf.ruled(ws.get(curve0), ws.get(curve1)))

    @add("translate", "Copy an object translated by a vector.",
         {"name": _NAME, "source": _NAME, "vector": _VEC}, ["name", "source", "vector"])
    def _(ws, name, source, vector):
        return ws.put(name, ws.get(source).transformed(nf.Transform.translation(vector)))

    @add("rotate", "Copy an object rotated about an axis. Angle in degrees.",
         {"name": _NAME, "source": _NAME, "axis_point": _VEC, "axis_dir": _VEC,
          "angle_deg": {"type": "number"}}, ["name", "source", "axis_point", "axis_dir", "angle_deg"])
    def _(ws, name, source, axis_point, axis_dir, angle_deg):
        t = nf.Transform.rotation(axis_point, axis_dir, _deg(angle_deg))
        return ws.put(name, ws.get(source).transformed(t))

    @add("mirror", "Copy an object mirrored at a plane (point + normal).",
         {"name": _NAME, "source": _NAME, "plane_point": _VEC, "plane_normal": _VEC},
         ["name", "source", "plane_point", "plane_normal"])
    def _(ws, name, source, plane_point, plane_normal):
        return ws.put(name, ws.get(source).transformed(nf.Transform.mirror(plane_point, plane_normal)))

    @add("polar_array", "Create count rotated copies about an axis; stored as name_0..name_{count-1}.",
         {"name": _NAME, "source": _NAME, "axis_point": _VEC, "axis_dir": _VEC,
          "count": {"type": "integer"}, "total_angle_deg": {"type": "number"}},
         ["name", "source", "axis_point", "axis_dir", "count", "total_angle_deg"])
    def _(ws, name, source, axis_point, axis_dir, count, total_angle_deg):
        ts = nf.polar_pattern(axis_point, axis_dir, int(count), _deg(total_angle_deg))
        src = ws.get(source)
        return "; ".join(ws.put(f"{name}_{i}", src.transformed(t)) for i, t in enumerate(ts))

    @add("linear_array", "Create count translated copies; stored as name_0..name_{count-1}.",
         {"name": _NAME, "source": _NAME, "step": _VEC, "count": {"type": "integer"}},
         ["name", "source", "step", "count"])
    def _(ws, name, source, step, count):
        src = ws.get(source)
        return "; ".join(ws.put(f"{name}_{i}", src.transformed(t))
                         for i, t in enumerate(nf.linear_pattern(step, int(count))))

    @add("list_objects", "List all objects in the workspace.", {}, [])
    def _(ws):
        return "\n".join(f"{k}: {v!r}" for k, v in ws.objects.items()) or "(empty)"

    @add("closest_point", "Closest point of an object to a given point, with its parameters "
         "and the distance.", {"name": _NAME, "point": _VEC}, ["name", "point"])
    def _(ws, name, point):
        r = nf.closest_point(ws.get(name), point)
        params = f"u={r.u:.6g}" if isinstance(r, nf.CurvePoint) else f"u={r.u:.6g}, v={r.v:.6g}"
        p = r.point
        return f"{params}, point=({p.x:.6g}, {p.y:.6g}, {p.z:.6g}), distance={r.distance:.6g}"

    @add("interpolate_points", "Create a smooth curve through the given points (in order).",
         {"name": _NAME, "points": {"type": "array", "items": _VEC, "minItems": 2},
          "degree": {"type": "integer", "description": "default 3"}}, ["name", "points"])
    def _(ws, name, points, degree=3):
        pts = np.asarray(points, dtype=float)
        return ws.put(name, nf.interpolate_curve(pts, min(int(degree), len(pts) - 1)))

    @add("loft", "Create a smooth surface through several section curves (in the given order).",
         {"name": _NAME, "sections": {"type": "array", "items": _NAME, "minItems": 2},
          "degree": {"type": "integer", "description": "degree across the sections, default 3"}},
         ["name", "sections"])
    def _(ws, name, sections, degree=3):
        return ws.put(name, nf.loft([ws.get(s) for s in sections], int(degree)))

    @add("sweep", "Sweep a profile curve along a path curve. The profile must be drawn at the "
         "start of the path (typically in the plane normal to the path). Optional twist in "
         "degrees over the whole path.",
         {"name": _NAME, "profile": _NAME, "path": _NAME, "twist_deg": {"type": "number"}},
         ["name", "profile", "path"])
    def _(ws, name, profile, path, twist_deg=0.0):
        opt = nf.SweepOptions(twist=_deg(twist_deg))
        return ws.put(name, nf.sweep(ws.get(profile), ws.get(path), opt))

    @add("coons", "Create a surface filling four boundary curves that form a closed loop "
         "(any orientation).",
         {"name": _NAME, "c0": _NAME, "c1": _NAME, "d0": _NAME, "d1": _NAME},
         ["name", "c0", "c1", "d0", "d1"])
    def _(ws, name, c0, c1, d0, d1):
        return ws.put(name, nf.coons(ws.get(c0), ws.get(c1), ws.get(d0), ws.get(d1)))

    @add("extend", "Extend a curve at its start or end (or a surface across a side umin/umax/"
         "vmin/vmax) by a length. mode 'natural' continues the shape smoothly, 'linear' "
         "continues straight along the tangent.",
         {"name": _NAME, "source": _NAME,
          "where": {"type": "string", "enum": ["start", "end", "umin", "umax", "vmin", "vmax"]},
          "length": {"type": "number"},
          "mode": {"type": "string", "enum": ["natural", "linear"]}},
         ["name", "source", "where", "length"])
    def _(ws, name, source, where, length, mode="natural"):
        o = ws.get(source)
        kind = nf.Extension.Linear if mode == "linear" else nf.Extension.Natural
        if isinstance(o, nf.NurbsCurve):
            end = {"start": nf.CurveEnd.Start, "end": nf.CurveEnd.End}[where]
            return ws.put(name, nf.extend(o, end, length, kind))
        side = {"umin": nf.SurfaceSide.UMin, "umax": nf.SurfaceSide.UMax,
                "vmin": nf.SurfaceSide.VMin, "vmax": nf.SurfaceSide.VMax}[where]
        return ws.put(name, nf.extend(o, side, length, kind))

    @add("measure", "Length of a curve or area of a surface.", {"name": _NAME}, ["name"])
    def _(ws, name):
        o = ws.get(name)
        if isinstance(o, nf.NurbsCurve):
            return f"length={nf.length(o):.10g}"
        return f"area={nf.area(o):.10g}"

    @add("evaluate", "Evaluate a curve at u or a surface at (u, v).",
         {"name": _NAME, "u": {"type": "number"}, "v": {"type": "number"}}, ["name", "u"])
    def _(ws, name, u, v=None):
        o = ws.get(name)
        return repr(o(u) if v is None else o(u, v))

    return T


TOOLS: list[Tool] = _build_tools()
TOOL_MAP: dict[str, Tool] = {t.name: t for t in TOOLS}


def call_tool(ws: Workspace, name: str, arguments: dict[str, Any]) -> str:
    """Execute a tool call; errors are returned as text so the LLM can correct itself."""
    tool = TOOL_MAP.get(name)
    if tool is None:
        return f"ERROR: unknown tool '{name}'"
    try:
        return tool.func(ws, **arguments)
    except Exception as exc:  # noqa: BLE001 - report every failure back to the model
        return f"ERROR: {type(exc).__name__}: {exc}"
