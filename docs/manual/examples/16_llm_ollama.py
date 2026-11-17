# requires: ollama
"""Natural-language geometry with a local model (needs `ollama serve` and a tool model)."""
from nurbforge.llm import GeometryAgent, OllamaClient

client = OllamaClient(host="http://localhost:11434", model="qwen2.5:14b")
print("models available:", client.list_models())
agent = GeometryAgent(client, verify_rounds=1)
print(agent.ask("Create a box of 20 x 10 x 5 mm at the origin and drill a hole of "
                "diameter 4 mm through its centre along z. Call the result 'plate'."))
for name, args, result in agent.log:
    print(f"  [{name}] {args} -> {result}")
print(agent.verification)
plate = agent.workspace.get("plate")       # an ordinary nurbforge Shell
print("volume:", plate.volume())
