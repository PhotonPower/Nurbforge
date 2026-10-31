# ADR 0004 – LLM access only through whitelisted tools

**Context.** Letting a model generate and execute Python would be flexible but unsafe and hard to
validate.

**Decision.** The Ollama agent can only call a fixed set of JSON-schema tools that map 1:1 to API
functions and act on a named `Workspace`. Errors are returned as text to the model.

**Consequences.** Safe for unattended use; capabilities grow by adding tools. Results are ordinary
Nurbforge objects that can be validated, meshed or exported.
