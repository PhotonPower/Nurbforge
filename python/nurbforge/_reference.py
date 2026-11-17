"""Compact API reference for language models and coding agents (``REFERENCE.md``)."""
from __future__ import annotations

from importlib import resources


def reference() -> str:
    """The compact API reference for language models and coding agents (Markdown): conventions,
    a worked example and every public function, class and LLM tool with signature and
    description. ``python -m nurbforge`` prints it."""
    return resources.files("nurbforge").joinpath("REFERENCE.md").read_text(encoding="utf-8")
