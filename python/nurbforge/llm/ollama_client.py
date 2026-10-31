"""Minimal Ollama HTTP client (stdlib only, no extra dependencies)."""
from __future__ import annotations

import json
import os
import urllib.error
import urllib.request
from typing import Any

DEFAULT_HOST = os.environ.get("OLLAMA_HOST", "http://localhost:11434")
DEFAULT_MODEL = os.environ.get("NURBFORGE_LLM_MODEL", "qwen2.5:14b")


class OllamaError(RuntimeError):
    pass


class OllamaClient:
    """Thin wrapper around the Ollama REST API (``/api/chat``, ``/api/tags``)."""

    def __init__(self, host: str = DEFAULT_HOST, model: str = DEFAULT_MODEL, timeout: float = 300.0):
        if not host.startswith("http"):
            host = "http://" + host
        self.host = host.rstrip("/")
        self.model = model
        self.timeout = timeout

    def _request(self, path: str, payload: dict[str, Any] | None = None) -> dict[str, Any]:
        data = json.dumps(payload).encode() if payload is not None else None
        req = urllib.request.Request(
            self.host + path, data=data, headers={"Content-Type": "application/json"},
            method="POST" if data is not None else "GET",
        )
        try:
            with urllib.request.urlopen(req, timeout=self.timeout) as resp:
                return json.loads(resp.read().decode())
        except urllib.error.URLError as exc:
            raise OllamaError(f"cannot reach Ollama at {self.host}: {exc}") from exc

    def list_models(self) -> list[str]:
        return [m["name"] for m in self._request("/api/tags").get("models", [])]

    def chat(self, messages: list[dict[str, Any]], tools: list[dict[str, Any]] | None = None,
             temperature: float = 0.0) -> dict[str, Any]:
        """Single non-streaming chat turn. Returns the assistant ``message`` dict."""
        payload: dict[str, Any] = {
            "model": self.model,
            "messages": messages,
            "stream": False,
            "options": {"temperature": temperature},
        }
        if tools:
            payload["tools"] = tools
        resp = self._request("/api/chat", payload)
        if "error" in resp:
            raise OllamaError(resp["error"])
        return resp["message"]
