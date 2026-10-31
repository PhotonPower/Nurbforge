# ADR 0002 – Dependency-free kernel

**Context.** Linear algebra libraries (Eigen) or full CAD kernels (OpenCASCADE) would speed up
development but complicate builds on AlmaLinux and Windows and embedding into solvers.

**Decision.** The core uses only the C++ standard library with small own `Vec3`/`Vec4`/`Transform`
types. Larger linear systems (interpolation, fitting) may later introduce Eigen as a *header-only,
private* dependency.

**Consequences.** Trivial builds and packaging; intersection and Boolean algorithms must be
implemented in-house, which is the core purpose of the project.
