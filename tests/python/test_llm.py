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
        self.seen_messages = list(messages)
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
        {"role": "assistant", "content": "VERIFIED: cylinder 'cyl' r=2, h=5."},
    ]
    agent = GeometryAgent(client=FakeClient(script))
    answer = agent.ask("Make a cylinder r=2, h=5")
    assert answer.startswith("VERIFIED") and "cyl" in answer
    assert "cyl: type=surface" in agent.verification and "bbox=" in agent.verification
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


def test_offset_tool():
    ws = Workspace()
    call_tool(ws, "make_circle", {"name": "c", "center": [0, 0, 0], "normal": [0, 0, 1], "radius": 2})
    out = call_tool(ws, "offset", {"name": "c2", "source": "c", "distance": -0.5})
    assert out.startswith("stored c2")
    assert abs(nf.norm(ws.get("c2")(0.3)) - 2.5) < 1e-7  # outwards for a ccw circle
    assert call_tool(ws, "offset", {"name": "bad", "source": "c", "distance": 2}).startswith("ERROR")
    call_tool(ws, "make_rectangle", {"name": "p", "origin": [0, 0, 0], "edge_u": [1, 0, 0],
                                     "edge_v": [0, 1, 0]})
    assert call_tool(ws, "offset", {"name": "p2", "source": "p", "distance": 0.25}).startswith("stored p2")
    assert abs(ws.get("p2")(0.5, 0.5).z - 0.25) < 1e-12


def test_fillet_tool():
    ws = Workspace()
    call_tool(ws, "make_line", {"name": "a", "p0": [-2, 0, 0], "p1": [0, 0, 0]})
    call_tool(ws, "make_line", {"name": "b", "p0": [0, 0, 0], "p1": [0, 2, 0]})
    out = call_tool(ws, "fillet", {"name": "f", "a": "a", "b": "b", "radius": 0.5})
    assert out.startswith("stored f") and "f_a" in out
    arc = ws.get("f")
    assert nf.distance(arc(arc.domain[0]), nf.Vec3(-0.5, 0, 0)) < 1e-12
    assert nf.distance(ws.get("f_b")(ws.get("f_b").domain[0]), nf.Vec3(0, 0.5, 0)) < 1e-12
    assert call_tool(ws, "fillet", {"name": "g", "a": "a", "b": "b", "radius": 9}).startswith("ERROR")
    call_tool(ws, "make_rectangle", {"name": "p", "origin": [-1, 0, 0], "edge_u": [5, 0, 0],
                                     "edge_v": [0, 3, 0]})
    call_tool(ws, "make_rectangle", {"name": "q", "origin": [0, 0, -1], "edge_u": [0, 3, 0],
                                     "edge_v": [0, 0, 5]})
    out = call_tool(ws, "fillet", {"name": "r", "a": "p", "b": "q", "radius": 0.5})
    assert out.startswith("stored r") and "r_contact1" in out
    c1 = ws.get("r_contact1")
    assert abs(c1(c1.domain[0]).z) < 1e-9
    assert call_tool(ws, "fillet", {"name": "x", "a": "a", "b": "p", "radius": 0.5}).startswith("ERROR")


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


def test_asphere_tool():
    ws = Workspace()
    out = call_tool(ws, "asphere", {"name": "lens", "radius": 50, "conic": -1, "clear_radius": 10})
    assert out.startswith("stored lens")
    p = ws.get("lens")(0.2, 1.0)
    assert abs(p.z - (p.x**2 + p.y**2) / 100.0) < 1e-12  # paraboloid


def test_intersect_tool():
    ws = Workspace()
    call_tool(ws, "make_line", {"name": "a", "p0": [0, 0, 0], "p1": [2, 2, 0]})
    call_tool(ws, "make_line", {"name": "b", "p0": [0, 2, 0], "p1": [2, 0, 0]})
    assert call_tool(ws, "intersect", {"a": "a", "b": "b"}) == "(1, 1, 0)"


def test_intersect_tool_surfaces():
    ws = Workspace()
    call_tool(ws, "make_rectangle", {"name": "p", "origin": [-2, -2, 0.5], "edge_u": [4, 0, 0],
                                     "edge_v": [0, 4, 0]})
    call_tool(ws, "make_circle", {"name": "c", "center": [0, 0, 0], "normal": [0, 0, 1], "radius": 1})
    call_tool(ws, "extrude", {"name": "cyl", "profile": "c", "direction": [0, 0, 1]})
    out = call_tool(ws, "intersect", {"a": "cyl", "b": "p"})
    assert "cut_1" in out and "length" in out


def test_trim_and_split_tools():
    ws = Workspace()
    call_tool(ws, "make_rectangle", {"name": "p", "origin": [0, 0, 0], "edge_u": [1, 0, 0], "edge_v": [0, 1, 0]})
    out = call_tool(ws, "trim_circle", {"name": "d", "source": "p", "center_u": 0.5, "center_v": 0.5,
                                        "radius": 0.5})
    assert out.startswith("stored d")
    assert call_tool(ws, "measure", {"name": "d"}).startswith("area=0.785398")
    call_tool(ws, "make_rectangle", {"name": "q", "origin": [0.5, -1, -1], "edge_u": [0, 3, 0],
                                     "edge_v": [0, 0, 2]})
    assert "s_2" in call_tool(ws, "split", {"name": "s", "surface": "p", "cutter": "q"})


def test_section_tool():
    ws = Workspace()
    call_tool(ws, "make_circle", {"name": "c", "center": [0, 0, 0], "normal": [0, 0, 1], "radius": 1})
    call_tool(ws, "extrude", {"name": "cyl", "profile": "c", "direction": [0, 0, 2]})
    out = call_tool(ws, "section", {"name": "cut", "surface": "cyl", "point": [0, 0, 1], "normal": [0, 0, 1]})
    assert "cut_1" in out


def test_sew_tool():
    ws = Workspace()
    corners = [((0, 0, 0), (1, 0, 0), (0, 1, 0)), ((0, 0, 0), (1, 0, 0), (0, 0, 1)),
               ((0, 1, 0), (0, 0, 1), (1, 0, 0)), ((0, 0, 0), (0, 1, 0), (0, 0, 1)),
               ((1, 0, 0), (0, 0, 1), (0, 1, 0)), ((0, 0, 1), (0, 1, 0), (1, 0, 0))]
    for i, (o, a, b) in enumerate(corners):
        call_tool(ws, "make_rectangle", {"name": f"f{i}", "origin": list(o), "edge_u": list(a),
                                         "edge_v": list(b)})
    out = call_tool(ws, "sew", {"name": "cube", "faces": [f"f{i}" for i in range(6)]})
    assert "closed shell, volume 1" in out


def test_boolean_tool():
    ws = Workspace()
    call_tool(ws, "make_arc", {"name": "m", "center": [0, 0, 0], "x_axis": [0, 0, -1], "y_axis": [1, 0, 0],
                               "radius": 1, "start_deg": 0, "end_deg": 180})
    call_tool(ws, "revolve", {"name": "s1", "profile": "m", "axis_point": [0, 0, 0], "axis_dir": [0, 0, 1],
                              "angle_deg": 360})
    call_tool(ws, "translate", {"name": "s2", "source": "s1", "vector": [1, 0, 0]})
    call_tool(ws, "sew", {"name": "a", "faces": ["s1"]})
    call_tool(ws, "sew", {"name": "b", "faces": ["s2"]})
    out = call_tool(ws, "boolean", {"name": "lens", "a": "a", "b": "b", "op": "intersection"})
    assert out.startswith("stored lens: closed shell, volume 1.30")


def test_inspect_and_check_tools():
    ws = Workspace()
    for name, lo, d1, d2 in [("b0", [0, 0, 0], [2, 0, 0], [0, 3, 0]), ("b1", [0, 0, 4], [2, 0, 0], [0, 3, 0]),
                             ("s0", [0, 0, 0], [2, 0, 0], [0, 0, 4]), ("s1", [0, 3, 0], [2, 0, 0], [0, 0, 4]),
                             ("s2", [0, 0, 0], [0, 3, 0], [0, 0, 4]), ("s3", [2, 0, 0], [0, 3, 0], [0, 0, 4])]:
        ws.put(name, nf.make_rectangle(lo, d1, d2))
    call_tool(ws, "sew", {"name": "box", "faces": ["b0", "b1", "s0", "s1", "s2", "s3"]})
    out = call_tool(ws, "inspect", {"name": "box"})
    assert "type=shell" in out and "closed=True" in out and "volume=24" in out and "ISSUES" not in out
    assert call_tool(ws, "check", {"name": "box", "property": "volume", "expected": 24}).startswith("PASS")
    assert call_tool(ws, "check", {"name": "box", "property": "bbox.max",
                                   "expected": [2, 3, 4]}).startswith("PASS")
    assert call_tool(ws, "check", {"name": "box", "property": "volume", "expected": 25}).startswith("FAIL")
    assert call_tool(ws, "check", {"name": "box", "property": "closed", "expected": True}).startswith("PASS")
    assert call_tool(ws, "check", {"name": "box", "property": "colour", "expected": 1}).startswith("ERROR")

    call_tool(ws, "sew", {"name": "lid", "faces": ["b0"]})
    assert "ISSUES: shell is open (4 free edges)" in call_tool(ws, "inspect", {"name": "lid"})


def test_verification_round_lets_the_model_fix_a_result():
    script = [
        # first attempt: a degenerate edge (p0 == p1)
        {"role": "assistant", "content": "", "tool_calls": [
            _call("make_line", name="edge", p0=[1, 0, 0], p1=[1, 0, 0])]},
        {"role": "assistant", "content": "Created 'edge'."},
        # verification: the report flags the issue, the model repairs the line
        {"role": "assistant", "content": "", "tool_calls": [
            _call("make_line", name="edge", p0=[0, 0, 0], p1=[1, 0, 0])]},
        {"role": "assistant", "content": "Fixed 'edge' (length 1)."},
        # second verification round: confirmed without tool calls
        {"role": "assistant", "content": "VERIFIED: 'edge' has length 1."},
    ]
    client = FakeClient(script)
    agent = GeometryAgent(client=client, verify_rounds=2)
    answer = agent.ask("Make a line of length 1 along x")
    assert answer.startswith("VERIFIED")
    assert not client.script
    first_report = [m["content"] for m in agent.messages if m["role"] == "user"][1]
    assert "ISSUES: degenerate curve (zero length)" in first_report
    assert "length=1" in agent.verification and "ISSUES" not in agent.verification
    assert abs(nf.length(agent.workspace.get("edge")) - 1.0) < 1e-12


def test_no_verification_without_new_objects():
    client = FakeClient([{"role": "assistant", "content": "Nothing to build."}])
    agent = GeometryAgent(client=client)
    assert agent.ask("hello") == "Nothing to build."
    assert agent.verification is None
