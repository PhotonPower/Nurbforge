# ADR 0001 – C++20 core with pybind11 Python API

**Context.** The kernel must be fast and embeddable into solver codes (C++/Fortran), yet scientists
mostly script in Python.

**Decision.** Implement all geometry in C++20 (`nf` namespace, static library). Expose it with
pybind11 as `nurbforge._core`; build Python wheels with scikit-build-core from the same CMake files.

**Consequences.** One source of truth for algorithms; Python adds only convenience, meshing and LLM
layers. Bindings must be kept in sync with the headers (enforced by Python tests).
nanobind remains an option should binary size or call overhead become relevant.
