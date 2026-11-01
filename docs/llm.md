# LLM interface (Ollama)

Nurbforge can build geometry from natural language using a local model served by
[Ollama](https://ollama.com). The interface uses Ollama's native **tool calling**.

## Safety model

The model **never executes code**. It can only call the whitelisted tools in
`nurbforge/llm/tools.py`, each with a JSON schema. Tools operate on a named object store
(`Workspace`). Invalid calls are returned to the model as `ERROR: ...` text so it can correct
itself; they never raise into the host program. This makes the interface usable in automated
pipelines.

## Usage

```bash
ollama serve &
ollama pull qwen2.5:14b          # any model with tool support (qwen2.5, llama3.1, mistral-nemo, ...)
nurbforge-llm --model qwen2.5:14b --out model.msh
```

```python
from nurbforge.llm import GeometryAgent, OllamaClient
agent = GeometryAgent(OllamaClient(host="http://localhost:11434", model="qwen2.5:14b"))
agent.ask("Revolve a semicircle of radius 5 about z into a sphere named ball")
ball = agent.workspace.get("ball")       # a regular NurbsSurface
```

Environment variables: `OLLAMA_HOST`, `NURBFORGE_LLM_MODEL`.

## Available tools

`make_line`, `make_circle`, `make_arc`, `make_rectangle`, `revolve`, `extrude`, `ruled`,
`translate`, `rotate`, `mirror`, `linear_array`, `polar_array`, `list_objects`, `evaluate`,
`closest_point`, `measure`, `interpolate_points`, `extend`, `loft`, `sweep`.
Angles are in **degrees** at this level (models handle degrees more reliably).

## Adding a tool

Add a function decorated with `@add(name, description, properties, required)` inside
`_build_tools()`. Keep descriptions short and precise; small models depend on them.
Add a test in `tests/python/test_llm.py` using the scripted `FakeClient`.
