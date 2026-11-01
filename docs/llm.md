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

## Verification

Language models are good at choosing operations and poor at noticing that a result is wrong –
a shell that did not close, a radius off by a factor of two, a part in the wrong place. Two
mechanisms let the model check its own work against measurements instead of guessing:

- **Tools** `inspect` (measured properties of an object: type, bounding box from dense samples,
  length / area / volume, closed, free and non-manifold edges, plus detected issues) and `check`
  (a measured property against an expected value, e.g. `{"name": "box", "property": "volume",
  "expected": 24}` → `PASS: volume = 24 (expected 24)`; properties like `bbox.max`, `closed`,
  `faces`; numbers within 1e-6 relative or a given absolute tolerance).
- **Verification rounds** (`GeometryAgent(verify_rounds=1)`, default; `nurbforge-llm
  --verify-rounds N`, 0 disables): after the model's final answer, every object created or
  changed for the request is measured and the report is sent back to the model. Issues are
  flagged explicitly – degenerate curves or surfaces, open shells with their number of free
  edges, non-manifold edges, inconsistent orientation, negative volume. The model then repairs
  the geometry with further tool calls or confirms (answer starting with `VERIFIED`); rounds
  stop as soon as it confirms without calling tools. The last report is kept in
  `agent.verification`.

```text
edge: type=curve, degree=1, ..., length=0, closed=True, ...  ISSUES: degenerate curve (zero length)
box: type=shell, faces=6, edges=12, free_edges=0, ..., closed=True, volume=24, bbox={min: [0, 0, 0], max: [2, 3, 4]}
```

The measurement code is in `nurbforge.llm.verify` (`inspect_object`, `issues`, `report`,
`check`) and can be used without a model.

## Available tools

`make_line`, `make_circle`, `make_arc`, `make_rectangle`, `revolve`, `extrude`, `ruled`,
`translate`, `rotate`, `mirror`, `linear_array`, `polar_array`, `list_objects`, `evaluate`,
`closest_point`, `measure`, `interpolate_points`, `extend`, `loft`, `sweep`, `coons`, `asphere`, `intersect`, `trim_circle`, `split`, `section`, `sew`, `boolean`, `inspect`, `check`,
`export_step` (native STEP writer; writes the file path given by the model).
Angles are in **degrees** at this level (models handle degrees more reliably).

## Adding a tool

Add a function decorated with `@add(name, description, properties, required)` inside
`_build_tools()`. Keep descriptions short and precise; small models depend on them.
Add a test in `tests/python/test_llm.py` using the scripted `FakeClient`.
