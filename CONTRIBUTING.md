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

## Documentation

The Markdown files in `docs/` are the source for both GitHub and the
[documentation site](https://photonpower.github.io/Nurbforge/) (Sphinx + MyST, furo theme). New
pages must be added to a `toctree` in `docs/index.md`. Build the site locally from the repository
root:

```bash
pip install -e ".[docs]"     # Sphinx extensions + compiled nurbforge._core (needed by autodoc)
doxygen docs/Doxyfile        # optional: C++ reference; writes build/docs/{html,xml}
sphinx-build -b html docs build/sphinx/html
```

Without Doxygen the build still succeeds; the C++ reference page then shows a placeholder.
Python docstrings and Doxygen comments in `include/` appear in the API reference, so keep them
current. The site is rebuilt and deployed by `.github/workflows/docs.yml` on every push to `main`.

## Packaging

Wheels are built with [cibuildwheel](https://cibuildwheel.pypa.io) (config in
`[tool.cibuildwheel]`, `pyproject.toml`) by `.github/workflows/wheels.yml` on `v*` tags, manual
dispatch and PRs touching `pyproject.toml`/`CMakeLists.txt`. To reproduce one locally (Windows,
no Docker needed): `pipx run cibuildwheel --only cp312-win_amd64`; Linux builds need Docker.

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
4. Document it in `docs/api/*.md` (docstrings/Doxygen comments feed the API reference), update
   `CHANGELOG.md`.

## Commit messages

Conventional style: `feat(construct): add loft surface`, `fix(basis): ...`, `docs: ...`.
