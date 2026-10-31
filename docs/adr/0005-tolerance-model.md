# ADR 0005 – Tolerance model

**Context.** Version 0.1 used one constant `kTolerance = 1e-10` for model-space lengths *and*
knot-value comparisons, plus ad-hoc literals (`1e-9`, `1e-12`) for angles. Intersections,
trimming and Booleans (roadmap 0.4/0.5) need clearly separated tolerances, because model-space
distances, parameter differences and angles have different units and scale differently.

**Decision.**

1. Three named kernel constants in `nurbforge/core/tolerance.hpp`:

   | Constant | Value | Unit | Used for |
   |---|---|---|---|
   | `nf::tol::kLength` | 1e-10 | model units | zero lengths, coincident points, parallel lines |
   | `nf::tol::kParametric` | 1e-10 | parameter units | equal knots, multiplicities, parameter snapping |
   | `nf::tol::kAngular` | 1e-9 | rad / dimensionless | angle limits, orthogonality of unit vectors |

   Weights are compared *relatively* (`|w − w₀| > kLength·|w₀|`), since a NURBS is invariant
   under uniform scaling of its weights.
2. No magic tolerance literals in `src/`; every comparison uses one of the constants above.
3. Iterative and approximating algorithms (point inversion, intersections, fitting, extrapolation)
   take a `const nf::Tolerance& tol = {}` argument whose defaults are the constants above. They
   document which fields they use and report the achieved accuracy where meaningful.
4. `nf::kTolerance` stays as an alias of `tol::kLength` for source compatibility.

**Consequences.** Validation behaviour of 0.1 is unchanged except that angular comparisons now
uniformly use 1e-9 (previously partly 1e-12) – harmless, as arcs remain exact. Models in very
large or very small units may need a per-operation `Tolerance`; a model-wide tolerance context
(stored on B-rep bodies, roadmap 0.5) is deferred until topology exists and will build on the
same struct.
