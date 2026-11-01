"""Tool-calling loop: natural language -> geometry operations."""
from __future__ import annotations

import json
from typing import Any

from .ollama_client import OllamaClient
from .tools import TOOLS, Workspace, call_tool
from .verify import report

SYSTEM_PROMPT = """You are the geometry assistant of Nurbforge, a NURBS CAD kernel.
Build geometry exclusively by calling the provided tools. Units are model units (e.g. mm).
Angles passed to tools are in degrees. Give every object a short, descriptive name and
reuse names of existing objects. Use 'inspect' and 'check' to verify dimensions, closedness
and volumes. If a dimension is only given indirectly (a target volume, area, length or
position), do not guess it: describe the construction as a recipe with unknowns and let 'solve'
find them. When finished, briefly summarise what was created."""

VERIFY_PROMPT = """Verification step. Measured properties of the objects created or changed for
this request:
{report}
Compare them with the request (dimensions, positions, closed solids, volumes). If anything is
wrong or listed under ISSUES, fix it with tool calls. Otherwise reply with a short final answer
starting with "VERIFIED"."""


class GeometryAgent:
    """Tool-calling agent.

    After the model's final answer, up to ``verify_rounds`` verification rounds follow: the
    objects created or changed for the request are measured (``verify.report``: bounding boxes,
    lengths, areas, volumes, closedness, detected issues such as open shells) and the report is
    sent back to the model, which either fixes the geometry with further tool calls or confirms.
    The rounds stop as soon as the model confirms without calling tools. The last report is kept
    in ``verification``."""

    def __init__(self, client: OllamaClient | None = None, workspace: Workspace | None = None,
                 max_steps: int = 20, verify_rounds: int = 1):
        self.client = client or OllamaClient()
        self.workspace = workspace or Workspace()
        self.max_steps = max_steps
        self.verify_rounds = verify_rounds
        self.messages: list[dict[str, Any]] = [{"role": "system", "content": SYSTEM_PROMPT}]
        self.log: list[tuple[str, dict[str, Any], str]] = []
        self.verification: str | None = None

    def ask(self, prompt: str) -> str:
        before = {n: id(o) for n, o in self.workspace.objects.items()}
        self.messages.append({"role": "user", "content": prompt})
        answer, _ = self._run()
        self.verification = None
        for _ in range(self.verify_rounds):
            changed = {n: o for n, o in self.workspace.objects.items() if before.get(n) != id(o)}
            if not changed:
                break
            self.verification = report(changed)
            self.messages.append({"role": "user",
                                  "content": VERIFY_PROMPT.format(report=self.verification)})
            answer, calls = self._run()
            if calls == 0:
                break
        return answer

    def _run(self) -> tuple[str, int]:
        """Chat until the model answers without tool calls; returns (answer, tool calls made)."""
        schemas = [t.schema() for t in TOOLS]
        n_calls = 0
        for _ in range(self.max_steps):
            msg = self.client.chat(self.messages, tools=schemas)
            self.messages.append(msg)
            calls = msg.get("tool_calls") or []
            if not calls:
                return msg.get("content", ""), n_calls
            for call in calls:
                fn = call["function"]
                args = fn.get("arguments") or {}
                if isinstance(args, str):
                    args = json.loads(args)
                result = call_tool(self.workspace, fn["name"], args)
                n_calls += 1
                self.log.append((fn["name"], args, result))
                self.messages.append({"role": "tool", "tool_name": fn["name"], "content": result})
        return "Stopped: maximum number of tool steps reached.", n_calls
