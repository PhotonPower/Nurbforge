"""Constraint solving for the LLM agent: dimensions derived from measured targets.

A request like "a cylinder of height 4 whose lateral area is 50" or "a lens of volume 2 cm³"
fixes a dimension only indirectly. Instead of guessing, the model describes the construction as
a *recipe* – a list of tool calls whose numeric arguments may refer to unknowns (``"$r"``) or
to arithmetic expressions over them (``"${2*r + 1}"``) – together with the unknowns (start
value and bounds) and *targets* on measured properties (``{"name": "cyl", "property": "area",
"value": 50}``, any property of ``verify.inspect_object``). ``solve`` runs the recipe in a
scratch copy of the workspace, measures the targets and adjusts the unknowns with a bounded
Levenberg–Marquardt iteration (finite-difference Jacobian). Only the final, converged recipe
is executed in the real workspace.

Expressions are evaluated by a small whitelist interpreter (numbers, unknowns, + - * / ** and
parentheses, pi, sqrt/sin/cos/tan/asin/acos/atan in degrees, abs, min, max) – the model never
executes Python code.
"""
from __future__ import annotations

import ast
import math
import re
from copy import copy
from typing import Any

from . import verify

_FUNCS = {
    "sqrt": math.sqrt, "abs": abs, "min": min, "max": max,
    "sin": lambda d: math.sin(math.radians(d)), "cos": lambda d: math.cos(math.radians(d)),
    "tan": lambda d: math.tan(math.radians(d)), "asin": lambda x: math.degrees(math.asin(x)),
    "acos": lambda x: math.degrees(math.acos(x)), "atan": lambda x: math.degrees(math.atan(x)),
}
_CONSTS = {"pi": math.pi}
_BINOPS = {ast.Add: lambda a, b: a + b, ast.Sub: lambda a, b: a - b, ast.Mult: lambda a, b: a * b,
           ast.Div: lambda a, b: a / b, ast.Pow: lambda a, b: a ** b}


def evaluate(expr: str, variables: dict[str, float]) -> float:
    """Evaluates an arithmetic expression over the variables (whitelist, no Python code)."""

    def ev(node: ast.AST) -> float:
        if isinstance(node, ast.Expression):
            return ev(node.body)
        if isinstance(node, ast.Constant) and isinstance(node.value, (int, float)) and not isinstance(
                node.value, bool):
            return float(node.value)
        if isinstance(node, ast.Name):
            if node.id in variables:
                return float(variables[node.id])
            if node.id in _CONSTS:
                return _CONSTS[node.id]
            raise ValueError(f"unknown name '{node.id}' in expression '{expr}'")
        if isinstance(node, ast.BinOp) and type(node.op) in _BINOPS:
            return _BINOPS[type(node.op)](ev(node.left), ev(node.right))
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.USub, ast.UAdd)):
            v = ev(node.operand)
            return -v if isinstance(node.op, ast.USub) else v
        if (isinstance(node, ast.Call) and isinstance(node.func, ast.Name) and node.func.id in _FUNCS
                and not node.keywords):
            return float(_FUNCS[node.func.id](*[ev(a) for a in node.args]))
        raise ValueError(f"unsupported expression '{expr}'")

    if len(expr) > 200:
        raise ValueError("expression too long")
    return ev(ast.parse(expr, mode="eval"))


_REF = re.compile(r"^\$(\w+)$|^\$\{(.+)\}$")


def substitute(value: Any, variables: dict[str, float]) -> Any:
    """Replaces "$name" and "${expression}" strings (also inside lists / dicts) by numbers."""
    if isinstance(value, str):
        m = _REF.match(value.strip())
        if m:
            return variables[m.group(1)] if m.group(1) else evaluate(m.group(2), variables)
        return value
    if isinstance(value, list):
        return [substitute(v, variables) for v in value]
    if isinstance(value, dict):
        return {k: substitute(v, variables) for k, v in value.items()}
    return value


def _run(ws: Any, recipe: list[dict[str, Any]], variables: dict[str, float], call_tool) -> None:
    for step in recipe:
        out = call_tool(ws, step["tool"], substitute(step.get("args", {}), variables))
        if out.startswith("ERROR"):
            raise RuntimeError(f"recipe step '{step['tool']}' failed: {out}")


def _residuals(ws: Any, recipe, variables, targets, call_tool) -> list[float]:
    scratch = copy(ws)
    scratch.objects = dict(ws.objects)
    _run(scratch, recipe, variables, call_tool)
    res = []
    for t in targets:
        value = verify._lookup(verify.inspect_object(scratch.get(t["name"])), t["property"])
        scale = max(abs(float(t["value"])), 1e-9)
        res.append((float(value) - float(t["value"])) / scale * float(t.get("weight", 1.0)))
    return res


def _solve_linear(A: list[list[float]], b: list[float]) -> list[float]:
    n = len(b)
    M = [row[:] + [b[i]] for i, row in enumerate(A)]
    for c in range(n):
        p = max(range(c, n), key=lambda r: abs(M[r][c]))
        if abs(M[p][c]) < 1e-300:
            raise ZeroDivisionError("singular system")
        M[c], M[p] = M[p], M[c]
        for r in range(c + 1, n):
            f = M[r][c] / M[c][c]
            for k in range(c, n + 1):
                M[r][k] -= f * M[c][k]
    x = [0.0] * n
    for r in range(n - 1, -1, -1):
        x[r] = (M[r][n] - sum(M[r][k] * x[k] for k in range(r + 1, n))) / M[r][r]
    return x


def solve(ws: Any, recipe: list[dict[str, Any]], unknowns: dict[str, dict[str, float]],
          targets: list[dict[str, Any]], call_tool, tolerance: float = 1e-8,
          max_iterations: int = 40) -> dict[str, Any]:
    """Bounded Levenberg–Marquardt on the relative target residuals; on success the recipe is
    executed in ``ws`` with the solution. Returns {"converged", "values", "residuals",
    "iterations"}."""
    names = list(unknowns)
    if not names or not targets:
        raise ValueError("solve needs at least one unknown and one target")
    lo = [float(unknowns[n].get("min", -math.inf)) for n in names]
    hi = [float(unknowns[n].get("max", math.inf)) for n in names]
    x = [min(max(float(unknowns[n]["start"]), lo[i]), hi[i]) for i, n in enumerate(names)]

    def F(v):
        return _residuals(ws, recipe, dict(zip(names, v)), targets, call_tool)

    def norm2(r):
        return sum(e * e for e in r)

    r = F(x)
    lam, it = 1e-3, 0
    while it < max_iterations and math.sqrt(norm2(r)) > tolerance:
        it += 1
        J = []  # forward differences, stepping inwards at the bounds
        for i in range(len(x)):
            h = 1e-6 * max(1.0, abs(x[i]))
            if x[i] + h > hi[i]:
                h = -h
            xp = x[:]
            xp[i] += h
            rp = F(xp)
            J.append([(a - b) / h for a, b in zip(rp, r)])
        # normal equations (J^T J + lam diag) dx = -J^T r ; J is stored by columns
        n = len(x)
        A = [[sum(J[i][k] * J[j][k] for k in range(len(r))) for j in range(n)] for i in range(n)]
        g = [sum(J[i][k] * r[k] for k in range(len(r))) for i in range(n)]
        improved = False
        for _ in range(12):
            Am = [[A[i][j] + (lam * max(A[i][i], 1e-12) if i == j else 0.0) for j in range(n)]
                  for i in range(n)]
            try:
                dx = _solve_linear(Am, [-gi for gi in g])
            except ZeroDivisionError:
                lam *= 10.0
                continue
            xn = [min(max(x[i] + dx[i], lo[i]), hi[i]) for i in range(n)]
            try:
                rn = F(xn)
            except Exception:  # noqa: BLE001 - an infeasible trial step is just rejected
                lam *= 10.0
                continue
            if norm2(rn) < norm2(r):
                x, r, lam, improved = xn, rn, max(lam / 10.0, 1e-12), True
                break
            lam *= 10.0
        if not improved:
            break
    converged = math.sqrt(norm2(r)) <= max(tolerance, 1e-6)
    values = dict(zip(names, x))
    if converged:
        _run(ws, recipe, values, call_tool)
    return {"converged": converged, "values": values, "residuals": r, "iterations": it}
