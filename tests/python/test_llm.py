"""LLM tool loop tested against a scripted fake Ollama client (no server needed)."""
import nurbforge as nf
from nurbforge.llm import GeometryAgent, Workspace
from nurbforge.llm.tools import TOOLS, call_tool


class FakeClient:
    def __init__(self, script):
        self.script = list(script)
        self.seen_tools = None

    def chat(self, messages, tools=None, temperature=0.0):
        self.seen_tools = tools
        return self.script.pop(0)


def _call(tool, **args):
    return {"function": {"name": tool, "arguments": args}}


def test_tool_schemas_are_valid():
    for t in TOOLS:
        s = t.schema()
        assert s["type"] == "function"
        assert s["function"]["parameters"]["type"] == "object"


def test_agent_builds_cylinder():
    script = [
        {"role": "assistant", "content": "", "tool_calls": [
            _call("make_line", name="gen", p0=[2, 0, 0], p1=[2, 0, 5]),
            _call("revolve", name="cyl", profile="gen", axis_point=[0, 0, 0],
                  axis_dir=[0, 0, 1], angle_deg=360),
        ]},
        {"role": "assistant", "content": "Created cylinder 'cyl'."},
    ]
    agent = GeometryAgent(client=FakeClient(script))
    answer = agent.ask("Make a cylinder r=2, h=5")
    assert "cyl" in answer
    cyl = agent.workspace.get("cyl")
    assert isinstance(cyl, nf.NurbsSurface)
    p = cyl(0.3, 1.0)
    assert abs((p.x**2 + p.y**2) ** 0.5 - 2.0) < 1e-12


def test_errors_are_reported_not_raised():
    ws = Workspace()
    out = call_tool(ws, "revolve", {"name": "x", "profile": "missing", "axis_point": [0, 0, 0],
                                    "axis_dir": [0, 0, 1], "angle_deg": 90})
    assert out.startswith("ERROR")
    assert call_tool(ws, "rm_rf", {}).startswith("ERROR")


def test_closest_point_tool():
    ws = Workspace()
    call_tool(ws, "make_circle", {"name": "c", "center": [0, 0, 0], "normal": [0, 0, 1], "radius": 1})
    out = call_tool(ws, "closest_point", {"name": "c", "point": [3, 0, 0]})
    assert "distance=2" in out


def test_measure_tool():
    ws = Workspace()
    call_tool(ws, "make_circle", {"name": "c", "center": [0, 0, 0], "normal": [0, 0, 1], "radius": 1})
    assert call_tool(ws, "measure", {"name": "c"}).startswith("length=6.28318530")


def test_interpolate_points_tool():
    ws = Workspace()
    out = call_tool(ws, "interpolate_points", {"name": "c", "points": [[0, 0, 0], [1, 1, 0], [2, 0, 0]]})
    assert out.startswith("stored c")
    assert ws.get("c").degree == 2  # degree limited by the number of points


def test_extend_tool():
    ws = Workspace()
    call_tool(ws, "make_line", {"name": "l", "p0": [0, 0, 0], "p1": [1, 0, 0]})
    out = call_tool(ws, "extend", {"name": "l2", "source": "l", "where": "end", "length": 1.5})
    assert out.startswith("stored l2")
    assert abs(ws.get("l2")(ws.get("l2").domain[1]).x - 2.5) < 1e-12


def test_loft_tool():
    ws = Workspace()
    for i, z in enumerate((0, 1, 2)):
        call_tool(ws, "make_circle", {"name": f"c{i}", "center": [0, 0, z], "normal": [0, 0, 1], "radius": 1})
    out = call_tool(ws, "loft", {"name": "tube", "sections": ["c0", "c1", "c2"]})
    assert out.startswith("stored tube")


def test_sweep_tool():
    ws = Workspace()
    call_tool(ws, "make_circle", {"name": "p", "center": [0, 0, 0], "normal": [0, 0, 1], "radius": 0.5})
    call_tool(ws, "make_line", {"name": "path", "p0": [0, 0, 0], "p1": [0, 0, 4]})
    out = call_tool(ws, "sweep", {"name": "pipe", "profile": "p", "path": "path", "twist_deg": 90})
    assert out.startswith("stored pipe")


def test_coons_tool():
    ws = Workspace()
    for n, a, b in (("c0", [0, 0, 0], [1, 0, 0]), ("c1", [0, 1, 0], [1, 1, 0]),
                    ("d0", [0, 0, 0], [0, 1, 0]), ("d1", [1, 0, 0], [1, 1, 0])):
        call_tool(ws, "make_line", {"name": n, "p0": a, "p1": b})
    out = call_tool(ws, "coons", {"name": "patch", "c0": "c0", "c1": "c1", "d0": "d0", "d1": "d1"})
    assert out.startswith("stored patch")
