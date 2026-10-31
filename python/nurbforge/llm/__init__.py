"""LLM interface (Ollama) for natural-language geometry construction.

Design principle: the model never executes arbitrary code. It can only call the
whitelisted tools in :mod:`nurbforge.llm.tools`, which operate on a named object
store (:class:`~nurbforge.llm.tools.Workspace`).
"""
from .agent import GeometryAgent  # noqa: F401
from .ollama_client import OllamaClient, OllamaError  # noqa: F401
from .tools import TOOLS, Workspace  # noqa: F401
