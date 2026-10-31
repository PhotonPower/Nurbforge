#pragma once
/// @file curve.hpp
/// Non-uniform rational B-spline curve in R^3.

#include <cstddef>
#include <vector>

#include "nurbforge/core/vec.hpp"
#include "nurbforge/nurbs/knot_vector.hpp"

namespace nf {

class Transform;

/// NURBS curve C(u) = sum N_i,p(u) w_i P_i / sum N_i,p(u) w_i.
///
/// Control points are stored in Cartesian form, weights separately.
/// Invariants (checked on construction): valid knot vector, weights > 0.
class NurbsCurve {
 public:
  /// Full constructor. Empty weights => non-rational (all weights 1).
  NurbsCurve(std::size_t degree, KnotVector knots, std::vector<Vec3> control_points,
             std::vector<double> weights = {});

  /// Non-rational B-spline with clamped uniform knots.
  static NurbsCurve from_points(std::size_t degree, std::vector<Vec3> control_points);

  [[nodiscard]] std::size_t degree() const noexcept { return degree_; }
  [[nodiscard]] const KnotVector& knots() const noexcept { return knots_; }
  [[nodiscard]] const std::vector<Vec3>& control_points() const noexcept { return ctrl_; }
  [[nodiscard]] const std::vector<double>& weights() const noexcept { return weights_; }
  [[nodiscard]] std::size_t num_control_points() const noexcept { return ctrl_.size(); }
  [[nodiscard]] bool is_rational() const noexcept;

  /// Parameter domain [u_p, u_{n+1}].
  [[nodiscard]] double domain_start() const noexcept { return knots_[degree_]; }
  [[nodiscard]] double domain_end() const noexcept { return knots_[ctrl_.size()]; }

  /// Point on the curve (A4.1).
  [[nodiscard]] Vec3 evaluate(double u) const;

  /// Derivatives C^(k)(u), k = 0..order (A3.2 + A4.2 for the rational case).
  [[nodiscard]] std::vector<Vec3> derivatives(double u, std::size_t order) const;

  /// Unit tangent vector at u.
  [[nodiscard]] Vec3 tangent(double u) const;

  /// Uniformly spaced samples in parameter space (n >= 2).
  [[nodiscard]] std::vector<Vec3> sample(std::size_t n) const;

  /// Returns a copy with knot u inserted r times (A5.1). Shape is unchanged.
  [[nodiscard]] NurbsCurve insert_knot(double u, std::size_t r = 1) const;

  /// Returns a copy with all knots in X inserted (A5.4). X may be unsorted and contain repeated
  /// values; every value must lie in the domain and no knot may exceed multiplicity p.
  /// Shape is unchanged.
  [[nodiscard]] NurbsCurve refine_knots(const std::vector<double>& X) const;

  /// Splits the curve into its Bézier segments (one per non-empty knot span of the domain).
  /// Segment k is a NURBS of the same degree with knots {a,..,a, b,..,b} on its original
  /// parameter interval [a, b], so parameters carry over unchanged.
  [[nodiscard]] std::vector<NurbsCurve> decompose_bezier() const;

  /// Returns the curve with reversed parametrisation (same geometry).
  [[nodiscard]] NurbsCurve reversed() const;

  /// Applies an affine transformation (exact: NURBS are affinely invariant).
  [[nodiscard]] NurbsCurve transformed(const Transform& t) const;

  /// Homogeneous control points (w*P, w).
  [[nodiscard]] std::vector<Vec4> homogeneous_control_points() const;

  /// Builds a curve from homogeneous control points.
  static NurbsCurve from_homogeneous(std::size_t degree, KnotVector knots,
                                     const std::vector<Vec4>& pw);

 private:
  std::size_t degree_;
  KnotVector knots_;
  std::vector<Vec3> ctrl_;
  std::vector<double> weights_;
};

}  // namespace nf
