"""Tool-calling loop: natural language -> geometry operations."""
from __future__ import annotations

import json
from typing import Any

from .ollama_client import OllamaClient
from .tools import TOOLS, Workspace, call_tool

SYSTEM_PROMPT = """You are the geometry assistant of Nurbforge, a NURBS CAD kernel.
Build geometry exclusively by calling the provided tools. Units are model units (e.g. mm).
Angles passed to tools are in degrees. Give every object a short, descriptive name and
reuse names of existing objects. When finished, briefly summarise what was created."""


class GeometryAgent:
    def __init__(self, client: OllamaClient | None = None, workspace: Workspace | None = None,
                 max_steps: int = 20):
        self.client = client or OllamaClient()
        self.workspace = workspace or Workspace()
        self.max_steps = max_steps
        self.messages: list[dict[str, Any]] = [{"role": "system", "content": SYSTEM_PROMPT}]
        self.log: list[tuple[str, dict[str, Any], str]] = []

    def ask(self, prompt: str) -> str:
        self.messages.append({"role": "user", "content": prompt})
        schemas = [t.schema() for t in TOOLS]
        for _ in range(self.max_steps):
            msg = self.client.chat(self.messages, tools=schemas)
            self.messages.append(msg)
            calls = msg.get("tool_calls") or []
            if not calls:
                return msg.get("content", "")
            for call in calls:
                fn = call["function"]
                args = fn.get("arguments") or {}
                if isinstance(args, str):
                    args = json.loads(args)
                result = call_tool(self.workspace, fn["name"], args)
                self.log.append((fn["name"], args, result))
                self.messages.append({"role": "tool", "tool_name": fn["name"], "content": result})
        return "Stopped: maximum number of tool steps reached."
