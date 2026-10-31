#include "nurbforge/nurbs/surface.hpp"

#include <algorithm>

#include "nurbforge/core/error.hpp"
#include "nurbforge/geometry/transform.hpp"
#include "nurbforge/nurbs/basis.hpp"
#include "nurbforge/nurbs/curve.hpp"

namespace nf {

NurbsSurface::NurbsSurface(std::size_t degree_u, std::size_t degree_v, KnotVector knots_u,
                           KnotVector knots_v, std::size_t num_u, std::size_t num_v,
                           std::vector<Vec3> control_points, std::vector<double> weights)
    : pu_(degree_u), pv_(degree_v), ku_(std::move(knots_u)), kv_(std::move(knots_v)),
      nu_(num_u), nv_(num_v), ctrl_(std::move(control_points)), weights_(std::move(weights)) {
  if (ctrl_.size() != nu_ * nv_) throw GeometryError("control net size must be num_u * num_v");
  if (weights_.empty()) weights_.assign(ctrl_.size(), 1.0);
  if (weights_.size() != ctrl_.size())
    throw GeometryError("number of weights must match number of control points");
  for (double w : weights_)
    if (!(w > 0.0)) throw GeometryError("weights must be strictly positive");
  validate_knots(ku_, pu_, nu_);
  validate_knots(kv_, pv_, nv_);
}

Vec3 NurbsSurface::evaluate(double u, double v) const {
  const std::size_t su = find_span(nu_ - 1, pu_, u, ku_);
  const std::size_t sv = find_span(nv_ - 1, pv_, v, kv_);
  const auto Nu = basis_functions(su, u, pu_, ku_);
  const auto Nv = basis_functions(sv, v, pv_, kv_);
  Vec4 sw;
  for (std::size_t a = 0; a <= pu_; ++a) {
    const std::size_t i = su - pu_ + a;
    Vec4 tmp;
    for (std::size_t b = 0; b <= pv_; ++b) {
      const std::size_t j = sv - pv_ + b;
      tmp += Nv[b] * Vec4::from_point(ctrl_[i * nv_ + j], weights_[i * nv_ + j]);
    }
    sw += Nu[a] * tmp;
  }
  return sw.project();
}

std::array<Vec3, 3> NurbsSurface::first_derivatives(double u, double v) const {
  const std::size_t su = find_span(nu_ - 1, pu_, u, ku_);
  const std::size_t sv = find_span(nv_ - 1, pv_, v, kv_);
  const auto Du = basis_function_derivatives(su, u, pu_, 1, ku_);
  const auto Dv = basis_function_derivatives(sv, v, pv_, 1, kv_);
  Vec3 A, Au, Av;
  double w = 0, wu = 0, wv = 0;
  for (std::size_t a = 0; a <= pu_; ++a) {
    for (std::size_t b = 0; b <= pv_; ++b) {
      const std::size_t idx = (su - pu_ + a) * nv_ + (sv - pv_ + b);
      const double wi = weights_[idx];
      const Vec3 pwi = wi * ctrl_[idx];
      A += Du[0][a] * Dv[0][b] * pwi;   w += Du[0][a] * Dv[0][b] * wi;
      Au += Du[1][a] * Dv[0][b] * pwi;  wu += Du[1][a] * Dv[0][b] * wi;
      Av += Du[0][a] * Dv[1][b] * pwi;  wv += Du[0][a] * Dv[1][b] * wi;
    }
  }
  const Vec3 S = A / w;
  return {S, (Au - wu * S) / w, (Av - wv * S) / w};
}

Vec3 NurbsSurface::normal(double u, double v) const {
  const auto d = first_derivatives(u, v);
  return normalized(cross(d[1], d[2]));
}

NurbsCurve NurbsSurface::iso_curve_u(double u) const {
  // Fix u: blend rows in homogeneous space, result runs along v.
  const std::size_t su = find_span(nu_ - 1, pu_, u, ku_);
  const auto Nu = basis_functions(su, u, pu_, ku_);
  std::vector<Vec4> pw(nv_);
  for (std::size_t j = 0; j < nv_; ++j)
    for (std::size_t a = 0; a <= pu_; ++a) {
      const std::size_t idx = (su - pu_ + a) * nv_ + j;
      pw[j] += Nu[a] * Vec4::from_point(ctrl_[idx], weights_[idx]);
    }
  return NurbsCurve::from_homogeneous(pv_, kv_, pw);
}

NurbsCurve NurbsSurface::iso_curve_v(double v) const { return transposed().iso_curve_u(v); }

std::vector<Vec3> NurbsSurface::sample(std::size_t nus, std::size_t nvs) const {
  if (nus < 2 || nvs < 2) throw GeometryError("sample() needs at least 2 samples per direction");
  const auto du = domain_u();
  const auto dv = domain_v();
  std::vector<Vec3> pts;
  pts.reserve(nus * nvs);
  for (std::size_t i = 0; i < nus; ++i) {
    const double u = du[0] + (du[1] - du[0]) * static_cast<double>(i) / static_cast<double>(nus - 1);
    for (std::size_t j = 0; j < nvs; ++j) {
      const double v = dv[0] + (dv[1] - dv[0]) * static_cast<double>(j) / static_cast<double>(nvs - 1);
      pts.push_back(evaluate(u, v));
    }
  }
  return pts;
}

NurbsSurface NurbsSurface::transformed(const Transform& t) const {
  std::vector<Vec3> c(ctrl_.size());
  std::transform(ctrl_.begin(), ctrl_.end(), c.begin(), [&](const Vec3& p) { return t.apply(p); });
  return {pu_, pv_, ku_, kv_, nu_, nv_, std::move(c), weights_};
}

NurbsSurface NurbsSurface::transposed() const {
  std::vector<Vec3> c(ctrl_.size());
  std::vector<double> w(ctrl_.size());
  for (std::size_t i = 0; i < nu_; ++i)
    for (std::size_t j = 0; j < nv_; ++j) {
      c[j * nu_ + i] = ctrl_[i * nv_ + j];
      w[j * nu_ + i] = weights_[i * nv_ + j];
    }
  return {pv_, pu_, kv_, ku_, nv_, nu_, std::move(c), std::move(w)};
}

}  // namespace nf
