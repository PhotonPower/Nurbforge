#pragma once
/// @file refine_detail.hpp
/// Internal knot refinement on homogeneous control points (shared by curves and surfaces).

#include <cstddef>
#include <utility>
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

/// Equivalent B-spline with clamped knot vector (first and last knot of multiplicity p+1) on the
/// same domain. Unchanged if already clamped. Knots closer than tol::kParametric are unified.
[[nodiscard]] RefinedPolygon clamp(std::size_t p, const KnotVector& U, const std::vector<Vec4>& Pw);

/// Degree elevation by t of a homogeneous B-spline (Piegl & Tiller A5.9). The input is clamped
/// first; the result is clamped, has degree p+t and every interior knot multiplicity raised by t.
[[nodiscard]] RefinedPolygon elevate_degree(std::size_t p, const KnotVector& U,
                                            const std::vector<Vec4>& Pw, std::size_t t);

/// Several control polygons sharing one degree and knot vector: a curve is one column, a surface
/// direction is the set of its columns (fixed index in the other direction).
struct Columns {
  KnotVector knots;
  std::vector<std::vector<Vec4>> cols;
};

/// Removes knot u up to `times` times from all columns simultaneously (A5.8), as long as every
/// column stays within the Cartesian deviation `tolerance` (converted with P&T Eq. 5.30).
/// Returns the number of removals done. Throws if u is not an interior knot.
std::size_t remove_knot(std::size_t p, Columns& net, double u, std::size_t times,
                        double tolerance);

/// Degree reduction by one of all columns, within the Cartesian deviation `tolerance`.
/// Bezier-wise reduction (P&T Eq. 5.41/5.42), then knot removal to restore the original
/// continuity; the result is verified rigorously by elevating it back and comparing control
/// polygons on a common knot vector. Throws GeometryError if the tolerance cannot be met.
[[nodiscard]] Columns reduce_degree(std::size_t p, const Columns& net, double tolerance);

/// Splits all columns at the interior parameter u (input is clamped first). Both parts are
/// clamped and keep their parameter intervals [a, u] and [u, b].
[[nodiscard]] std::pair<Columns, Columns> split(std::size_t p, const Columns& net, double u);

/// Clamps all columns (common knot vector).
[[nodiscard]] Columns clamp(std::size_t p, const Columns& net);

/// Start indices k of the non-empty spans [U[k], U[k+1]) inside the domain. After Bézier
/// refinement the segment on span k has control points k-p .. k.
[[nodiscard]] std::vector<std::size_t> nonempty_spans(std::size_t p, const KnotVector& U,
                                                      std::size_t n_ctrl);

}  // namespace nf::detail
