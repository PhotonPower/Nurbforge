#pragma once
/// @file knot_vector.hpp
/// Knot vector utilities (validation, span search, generation).

#include <cstddef>
#include <utility>
#include <vector>

namespace nf {

using KnotVector = std::vector<double>;

/// Validates a knot vector for a B-spline of given degree with n_ctrl control points.
/// Checks size (n_ctrl + degree + 1), monotonicity and interior multiplicity <= degree.
/// Throws GeometryError on failure.
void validate_knots(const KnotVector& knots, std::size_t degree, std::size_t n_ctrl);

/// Clamped (open) uniform knot vector on [0, 1].
[[nodiscard]] KnotVector clamped_uniform_knots(std::size_t degree, std::size_t n_ctrl);

/// Index of the knot span containing u (Piegl & Tiller, A2.1).
/// n = n_ctrl - 1. Returns n for u == knots[n+1].
[[nodiscard]] std::size_t find_span(std::size_t n, std::size_t degree, double u,
                                    const KnotVector& knots);

/// Multiplicity of value u inside the knot vector (within tolerance).
[[nodiscard]] std::size_t knot_multiplicity(const KnotVector& knots, double u);

/// Distinct knot values and their multiplicities (format used e.g. by gmsh / OCC).
[[nodiscard]] std::pair<std::vector<double>, std::vector<int>> unique_knots(const KnotVector& knots);

}  // namespace nf
