"""Geometry verification for the LLM agent: measured properties and sanity issues of objects.

``inspect_object`` measures an object of the workspace (curve, surface, trimmed face, shell);
``issues`` lists problems a model cannot see from tool messages alone – an open shell, an
inverted solid, degenerate geometry. ``report`` formats both for the model, and ``check``
compares a measured property with an expected value (used by the ``check`` tool).
"""
from __future__ import annotations

import math
from typing import Any

from .. import _core as nf

_REL = 1e-8  # relative accuracy of measured areas / volumes / lengths


def _bbox(points: list[Any]) -> dict[str, list[float]]:
    xs = [p.x for p in points]
    ys = [p.y for p in points]
    zs = [p.z for p in points]
    return {"min": [min(xs), min(ys), min(zs)], "max": [max(xs), max(ys), max(zs)]}


def _vec(p: Any) -> list[float]:
    return [p.x, p.y, p.z]


def _surface_points(s: Any, n: int = 16, keep=None) -> list[Any]:
    (u0, u1), (v0, v1) = s.domain_u, s.domain_v
    pts = []
    for i in range(n + 1):
        for j in range(n + 1):
            u = u0 + (u1 - u0) * i / n
            v = v0 + (v1 - v0) * j / n
            if keep is None or keep(u, v):
                pts.append(s(u, v))
    return pts


def inspect_object(obj: Any) -> dict[str, Any]:
    """Measured properties of a workspace object (bounding boxes from dense samples)."""
    if isinstance(obj, nf.NurbsCurve):
        a, b = obj.domain
        p0, p1 = obj(a), obj(b)
        pts = [obj(a + (b - a) * k / 64) for k in range(65)]
        return {"type": "curve", "degree": obj.degree, "control_points": len(obj.control_points),
                "rational": obj.is_rational, "length": nf.length(obj), "start": _vec(p0),
                "end": _vec(p1), "closed": nf.distance(p0, p1) <= 1e-9 * max(1.0, nf.length(obj)),
                "bbox": _bbox(pts)}
    if isinstance(obj, nf.NurbsSurface):
        return {"type": "surface", "degree": [obj.degree_u, obj.degree_v],
                "control_points": [obj.num_u, obj.num_v], "area": nf.area(obj),
                "bbox": _bbox(_surface_points(obj))}
    if isinstance(obj, nf.TrimmedSurface):
        pts = _surface_points(obj.surface, 24, obj.contains)
        out = {"type": "trimmed_surface", "holes": len(obj.holes), "area": obj.area(_REL)}
        if pts:
            out["bbox"] = _bbox(pts)
        return out
    if isinstance(obj, nf.Shell):
        pts = []
        for f in obj.faces:
            pts += _surface_points(f.surface, 12, f.contains)
        out = {"type": "shell", "faces": len(obj.faces), "edges": len(obj.edges),
               "free_edges": len(obj.free_edges()), "non_manifold_edges": len(obj.non_manifold_edges()),
               "closed": obj.is_closed, "consistently_oriented": obj.is_consistently_oriented,
               "area": obj.area(_REL)}
        if obj.is_closed and obj.is_consistently_oriented:
            out["volume"] = obj.volume(_REL)
        if pts:
            out["bbox"] = _bbox(pts)
        return out
    return {"type": type(obj).__name__}


def issues(props: dict[str, Any]) -> list[str]:
    """Sanity problems of an object, from its measured properties."""
    out = []
    t = props.get("type")
    if t == "curve" and props["length"] <= 1e-12:
        out.append("degenerate curve (zero length)")
    if t in ("surface", "trimmed_surface") and props["area"] <= 1e-12:
        out.append("degenerate surface (zero area)")
    if t == "shell":
        if props["faces"] == 0:
            out.append("empty shell (no faces)")
        elif not props["closed"]:
            out.append(f"shell is open ({props['free_edges']} free edges) – not a solid")
        if props["non_manifold_edges"]:
            out.append(f"{props['non_manifold_edges']} non-manifold edges")
        if props["closed"] and not props["consistently_oriented"]:
            out.append("faces are not consistently oriented")
        if props.get("volume", 1.0) < 0.0:
            out.append("negative volume (shell oriented inwards)")
    return out


def _fmt(v: Any) -> str:
    if isinstance(v, float):
        return f"{v:.10g}"
    if isinstance(v, list):
        return "[" + ", ".join(_fmt(x) for x in v) + "]"
    if isinstance(v, dict):
        return "{" + ", ".join(f"{k}: {_fmt(x)}" for k, x in v.items()) + "}"
    return str(v)


def describe(name: str, obj: Any) -> str:
    """One line per object: measured properties and issues."""
    try:
        props = inspect_object(obj)
    except Exception as exc:  # noqa: BLE001 - a failing measurement is itself a finding
        return f"{name}: could not be measured ({type(exc).__name__}: {exc})"
    line = f"{name}: " + ", ".join(f"{k}={_fmt(v)}" for k, v in props.items())
    found = issues(props)
    if found:
        line += "  ISSUES: " + "; ".join(found)
    return line


def report(objects: dict[str, Any]) -> str:
    """Verification report for several named objects."""
    return "\n".join(describe(n, o) for n, o in objects.items())


def _lookup(props: dict[str, Any], prop: str) -> Any:
    v: Any = props
    for part in prop.replace("[", ".").replace("]", "").split("."):
        if isinstance(v, dict):
            if part not in v:
                raise KeyError(f"no property '{prop}'. Available: {sorted(props)}")
            v = v[part]
        else:
            v = v[int(part)]
    return v


def check(obj: Any, prop: str, expected: Any, tolerance: float | None = None) -> str:
    """Compare a measured property (e.g. 'volume', 'bbox.max', 'closed') with an expectation.

    Numbers (and number lists) pass within ``tolerance`` (default: 1e-6 relative, at least
    1e-9 absolute); other values must be equal."""
    props = inspect_object(obj)
    value = _lookup(props, prop)

    def close(a: float, b: float) -> bool:
        tol = tolerance if tolerance is not None else max(1e-9, 1e-6 * abs(b))
        return math.isfinite(a) and abs(a - b) <= tol

    if isinstance(value, (bool, str)) or isinstance(expected, (bool, str)):
        ok = value == expected
    elif isinstance(value, list):
        ok = isinstance(expected, list) and len(value) == len(expected) and all(
            close(float(a), float(b)) for a, b in zip(value, expected))
    else:
        ok = close(float(value), float(expected))
    return f"{'PASS' if ok else 'FAIL'}: {prop} = {_fmt(value)} (expected {_fmt(expected)})"
