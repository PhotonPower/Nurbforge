# Security policy

Please report vulnerabilities privately via GitHub's *Report a vulnerability* (Security tab)
rather than in public issues.

Notes for users:
- The LLM interface only executes whitelisted geometry tools (see `docs/llm.md`); it does not run
  model-generated code.
- Never commit access tokens. If a token is exposed, revoke it immediately in GitHub settings.
