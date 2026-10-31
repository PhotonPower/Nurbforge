"""Natural-language geometry via a local Ollama model.

Prerequisites:  ollama serve   and   ollama pull qwen2.5:14b   (any tool-capable model)
"""
from nurbforge.llm import GeometryAgent, OllamaClient

agent = GeometryAgent(OllamaClient(model="qwen2.5:14b"))
print(agent.ask("Create a cylinder with radius 10 and height 30 along z, "
                "then a polar array of 6 copies of a circle of radius 2 at x=20 in the xy-plane."))
for name, args, result in agent.log:
    print(f"[{name}] {args} -> {result}")
