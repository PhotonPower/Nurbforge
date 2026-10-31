#pragma once
/// @file tolerance.hpp
/// Named tolerances of the kernel (see docs/adr/0005-tolerance-model.md).
///
/// Three independent quantities are distinguished because they have different units:
///  - length:     distances in model space (model units)
///  - parametric: differences of knot / parameter values (parameter units)
///  - angular:    angles in radians and dot products of unit vectors (dimensionless)
///
/// The constants are the kernel defaults used for input validation and exact comparisons.
/// Iterative algorithms (point inversion, intersections, fitting) take a `Tolerance` argument
/// so callers can loosen or tighten them per operation.

namespace nf {

namespace tol {
inline constexpr double kLength = 1e-10;      ///< model units
inline constexpr double kParametric = 1e-10;  ///< parameter units
inline constexpr double kAngular = 1e-9;      ///< radians / dimensionless
}  // namespace tol

/// Per-operation tolerance set. Defaults equal the kernel constants in `nf::tol`.
struct Tolerance {
  double length{tol::kLength};
  double parametric{tol::kParametric};
  double angular{tol::kAngular};
};

}  // namespace nf
