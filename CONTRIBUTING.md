# Contributing to Nurbforge

## Workflow

1. Open an issue describing the feature or bug.
2. Branch from `main`: `feature/<topic>` or `fix/<topic>`.
3. Add tests (C++: `tests/cpp`, Python: `tests/python`) – every geometric algorithm needs at least
   one test against an exact reference (e.g. points of a revolved circle lie on the sphere).
4. Run locally:
   ```bash
   cmake --preset debug && cmake --build --preset debug && ctest --preset debug
   pip install -e ".[dev]" && pytest
   ```
5. Open a pull request; CI must pass on AlmaLinux 10 and Windows 11.

## Code style

- C++20, formatted with `clang-format` (config in repo); `clang-tidy` config provided.
- Namespace `nf`; types `PascalCase`, functions/variables `snake_case`, members end in `_`.
- Geometry objects are immutable values; operations return new objects.
- Reference *The NURBS Book* algorithm numbers in comments where applicable.
- Python: `ruff` (line length 100), type hints, docstrings.

## Adding a feature end-to-end

1. Header in `include/nurbforge/<area>/`, implementation in `src/<area>/`, add to `CMakeLists.txt`.
2. Bind it in `python/bindings/module.cpp` and re-export in `python/nurbforge/__init__.py`.
3. Optional: expose to the LLM in `python/nurbforge/llm/tools.py`.
4. Document it in `docs/api/*.md`, update `CHANGELOG.md`.

## Commit messages

Conventional style: `feat(construct): add loft surface`, `fix(basis): ...`, `docs: ...`.
