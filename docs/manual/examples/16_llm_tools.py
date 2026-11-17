"""The LLM tool layer without a model: workspace, tool calls, verification, solving."""
import json

from nurbforge.llm import TOOLS, Workspace, verify
from nurbforge.llm.tools import call_tool

print(f"{len(TOOLS)} tools, e.g.", ", ".join(t.name for t in TOOLS[:6]), "...")
schema = TOOLS[1].schema()["function"]
print(schema["name"], "-", schema["description"])
print("  parameters:", json.dumps(schema["parameters"]["properties"]["radius"]),
      "required:", schema["parameters"]["required"])

ws = Workspace()
# every call returns a text message, exactly what the model would see
print(call_tool(ws, "solid", {"name": "box", "kind": "box", "point": [0, 0, 0],
                              "size": [2, 3, 4]}))
print(call_tool(ws, "solid", {"name": "rod", "kind": "cylinder", "point": [1, 1.5, -1],
                              "axis": [0, 0, 1], "radius": 0.5, "height": 6}))
print(call_tool(ws, "boolean", {"name": "part", "a": "box", "b": "rod", "op": "difference"}))
print(call_tool(ws, "make_arc", {"name": "a", "center": [0, 0, 0], "x_axis": [1, 0, 0],
                                 "y_axis": [0, 1, 0], "radius": 2, "start_deg": 0,
                                 "end_deg": 90}))
print(call_tool(ws, "measure", {"name": "a"}))

# errors never raise: they come back as text, so the model can correct itself
print(call_tool(ws, "boolean", {"name": "x", "a": "box", "b": "nothing", "op": "union"}))
print(call_tool(ws, "no_such_tool", {}))

# verification: measured properties and detected issues
print(call_tool(ws, "inspect", {"name": "part"}))
print(call_tool(ws, "check", {"name": "box", "property": "volume", "expected": 24}))
print(call_tool(ws, "check", {"name": "box", "property": "bbox.max", "expected": [2, 3, 5]}))
print(call_tool(ws, "make_rectangle", {"name": "plate", "origin": [0, 0, 0],
                                       "edge_u": [1, 0, 0], "edge_v": [0, 1, 0]}))
print(call_tool(ws, "sew", {"name": "lid", "faces": ["plate"]}))
print(verify.report({"lid": ws.get("lid"), "a": ws.get("a")}))

# constraint solving: a cylinder of height 4 with a lateral area of 50
recipe = [
    {"tool": "make_line", "args": {"name": "gen", "p0": ["$r", 0, 0], "p1": ["$r", 0, "$h"]}},
    {"tool": "revolve", "args": {"name": "cyl", "profile": "gen", "axis_point": [0, 0, 0],
                                 "axis_dir": [0, 0, 1], "angle_deg": 360}},
]
print(call_tool(ws, "solve", {
    "recipe": recipe,
    "unknowns": {"r": {"start": 1, "min": 0.01, "max": 10},
                 "h": {"start": 1, "min": 0.01, "max": 10}},
    "targets": [{"name": "cyl", "property": "area", "value": 50},
                {"name": "cyl", "property": "bbox.max[2]", "value": 4}]}))
print(call_tool(ws, "check", {"name": "cyl", "property": "area", "expected": 50}))
print("workspace:", sorted(ws.objects))
