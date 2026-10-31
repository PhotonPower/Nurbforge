#pragma once
/// @file surface.hpp
/// Tensor-product NURBS surface in R^3.

#include <array>
#include <cstddef>
#include <vector>

#include "nurbforge/core/vec.hpp"
#include "nurbforge/nurbs/knot_vector.hpp"

namespace nf {

class Transform;
class NurbsCurve;

/// NURBS surface S(u,v). Control net is stored row-major: index(i, j) = i * num_v + j,
/// where i runs along u (0..num_u-1) and j along v (0..num_v-1).
class NurbsSurface {
 public:
  NurbsSurface(std::size_t degree_u, std::size_t degree_v, KnotVector knots_u, KnotVector knots_v,
               std::size_t num_u, std::size_t num_v, std::vector<Vec3> control_points,
               std::vector<double> weights = {});

  [[nodiscard]] std::size_t degree_u() const noexcept { return pu_; }
  [[nodiscard]] std::size_t degree_v() const noexcept { return pv_; }
  [[nodiscard]] const KnotVector& knots_u() const noexcept { return ku_; }
  [[nodiscard]] const KnotVector& knots_v() const noexcept { return kv_; }
  [[nodiscard]] std::size_t num_u() const noexcept { return nu_; }
  [[nodiscard]] std::size_t num_v() const noexcept { return nv_; }
  [[nodiscard]] const std::vector<Vec3>& control_points() const noexcept { return ctrl_; }
  [[nodiscard]] const std::vector<double>& weights() const noexcept { return weights_; }
  [[nodiscard]] const Vec3& control_point(std::size_t i, std::size_t j) const { return ctrl_[i * nv_ + j]; }
  [[nodiscard]] double weight(std::size_t i, std::size_t j) const { return weights_[i * nv_ + j]; }

  [[nodiscard]] std::array<double, 2> domain_u() const noexcept { return {ku_[pu_], ku_[nu_]}; }
  [[nodiscard]] std::array<double, 2> domain_v() const noexcept { return {kv_[pv_], kv_[nv_]}; }

  /// Point on the surface (A4.3).
  [[nodiscard]] Vec3 evaluate(double u, double v) const;

  /// {S, S_u, S_v} – first order partial derivatives (rational).
  [[nodiscard]] std::array<Vec3, 3> first_derivatives(double u, double v) const;

  /// Unit normal (S_u x S_v). Throws at degenerate points (e.g. poles).
  [[nodiscard]] Vec3 normal(double u, double v) const;

  /// Iso-parametric curve at fixed u (direction v) or fixed v (direction u).
  [[nodiscard]] NurbsCurve iso_curve_u(double u) const;
  [[nodiscard]] NurbsCurve iso_curve_v(double v) const;

  /// Structured grid of points (nu_samples x nv_samples), row-major.
  [[nodiscard]] std::vector<Vec3> sample(std::size_t nu_samples, std::size_t nv_samples) const;

  /// Inserts all knots X into the u (resp. v) knot vector (A5.5). Shape is unchanged.
  /// Same rules as NurbsCurve::refine_knots.
  [[nodiscard]] NurbsSurface refine_knots_u(const std::vector<double>& X) const;
  [[nodiscard]] NurbsSurface refine_knots_v(const std::vector<double>& X) const;

  /// Splits the surface into Bézier patches, one per pair of non-empty knot spans. Returned in
  /// row-major order (u span outer, v span inner); each patch keeps its original parameter
  /// rectangle.
  [[nodiscard]] std::vector<NurbsSurface> decompose_bezier() const;

  [[nodiscard]] NurbsSurface transformed(const Transform& t) const;

  /// Swaps the u and v directions (transposes the control net).
  [[nodiscard]] NurbsSurface transposed() const;

 private:
  std::size_t pu_, pv_;
  KnotVector ku_, kv_;
  std::size_t nu_, nv_;
  std::vector<Vec3> ctrl_;
  std::vector<double> weights_;
};

}  // namespace nf
