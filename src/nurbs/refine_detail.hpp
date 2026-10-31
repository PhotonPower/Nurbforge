#pragma once
/// @file refine_detail.hpp
/// Internal knot refinement on homogeneous control points (shared by curves and surfaces).

#include <cstddef>
#include <vector>

#include "nurbforge/core/vec.hpp"
#include "nurbforge/nurbs/knot_vector.hpp"

namespace nf::detail {

/// Result of refining a homogeneous B-spline.
struct RefinedPolygon {
  KnotVector knots;
  std::vector<Vec4> points;
};

/// Inserts all knots X (any order) into the homogeneous B-spline (p, U, Pw) – Piegl & Tiller
/// A5.4. Validates that every x lies in the parameter domain and that no knot value would get a
/// multiplicity above p. Throws GeometryError otherwise.
[[nodiscard]] RefinedPolygon refine_knot_vector(std::size_t p, const KnotVector& U,
                                                const std::vector<Vec4>& Pw,
                                                std::vector<double> X);

/// Knots to insert so that every distinct knot value in the parameter domain has multiplicity
/// >= p, i.e. the B-spline splits into Bézier segments (sorted ascending).
[[nodiscard]] std::vector<double> bezier_refinement_knots(std::size_t p, const KnotVector& U,
                                                          std::size_t n_ctrl);

/// Start indices k of the non-empty spans [U[k], U[k+1]) inside the domain. After Bézier
/// refinement the segment on span k has control points k-p .. k.
[[nodiscard]] std::vector<std::size_t> nonempty_spans(std::size_t p, const KnotVector& U,
                                                      std::size_t n_ctrl);

}  // namespace nf::detail
