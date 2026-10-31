#include "nurbforge/nurbs/surface.hpp"

#include <algorithm>
#include <string>

#include "nurbforge/core/error.hpp"
#include "nurbforge/geometry/transform.hpp"
#include "nurbforge/nurbs/basis.hpp"
#include "nurbforge/nurbs/curve.hpp"
#include "refine_detail.hpp"

namespace nf {

namespace {

// Columns of the u direction: column j holds the points (i, j), i = 0..num_u-1.
detail::Columns columns_u(const NurbsSurface& s) {
  detail::Columns net{s.knots_u(), {}};
  for (std::size_t j = 0; j < s.num_v(); ++j) {
    std::vector<Vec4> col(s.num_u());
    for (std::size_t i = 0; i < s.num_u(); ++i)
      col[i] = Vec4::from_point(s.control_point(i, j), s.weight(i, j));
    net.cols.push_back(std::move(col));
  }
  return net;
}

NurbsSurface from_columns_u(const NurbsSurface& s, std::size_t degree_u, const detail::Columns& net) {
  const std::size_t nu = net.cols.front().size(), nv = net.cols.size();
  std::vector<Vec3> c(nu * nv);
  std::vector<double> w(nu * nv);
  for (std::size_t j = 0; j < nv; ++j)
    for (std::size_t i = 0; i < nu; ++i) {
      c[i * nv + j] = net.cols[j][i].project();
      w[i * nv + j] = net.cols[j][i].w;
    }
  return {degree_u, s.degree_v(), net.knots, s.knots_v(), nu, nv, std::move(c), std::move(w)};
}

}  // namespace

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

NurbsSurface NurbsSurface::refine_knots_u(const std::vector<double>& X) const {
  if (X.empty()) return *this;
  KnotVector new_ku;
  std::size_t new_nu = 0;
  std::vector<Vec3> c;
  std::vector<double> w;
  // refine every column (fixed j, i along u) – A5.5 applied column by column
  for (std::size_t j = 0; j < nv_; ++j) {
    std::vector<Vec4> col(nu_);
    for (std::size_t i = 0; i < nu_; ++i)
      col[i] = Vec4::from_point(ctrl_[i * nv_ + j], weights_[i * nv_ + j]);
    auto r = detail::refine_knot_vector(pu_, ku_, col, X);
    if (j == 0) {
      new_ku = std::move(r.knots);
      new_nu = r.points.size();
      c.resize(new_nu * nv_);
      w.resize(new_nu * nv_);
    }
    for (std::size_t i = 0; i < new_nu; ++i) {
      c[i * nv_ + j] = r.points[i].project();
      w[i * nv_ + j] = r.points[i].w;
    }
  }
  return {pu_, pv_, std::move(new_ku), kv_, new_nu, nv_, std::move(c), std::move(w)};
}

NurbsSurface NurbsSurface::elevate_degree_u(std::size_t t) const {
  if (t == 0) return *this;
  KnotVector new_ku;
  std::size_t new_nu = 0;
  std::vector<Vec3> c;
  std::vector<double> w;
  for (std::size_t j = 0; j < nv_; ++j) {  // every column (fixed j) has the same knots
    std::vector<Vec4> col(nu_);
    for (std::size_t i = 0; i < nu_; ++i)
      col[i] = Vec4::from_point(ctrl_[i * nv_ + j], weights_[i * nv_ + j]);
    auto r = detail::elevate_degree(pu_, ku_, col, t);
    if (j == 0) {
      new_ku = std::move(r.knots);
      new_nu = r.points.size();
      c.resize(new_nu * nv_);
      w.resize(new_nu * nv_);
    }
    for (std::size_t i = 0; i < new_nu; ++i) {
      c[i * nv_ + j] = r.points[i].project();
      w[i * nv_ + j] = r.points[i].w;
    }
  }
  return {pu_ + t, pv_, std::move(new_ku), kv_, new_nu, nv_, std::move(c), std::move(w)};
}

NurbsSurface NurbsSurface::elevate_degree_v(std::size_t t) const {
  if (t == 0) return *this;
  return transposed().elevate_degree_u(t).transposed();
}

NurbsSurface NurbsSurface::clamped() const {
  const bool cu = ku_.front() == ku_[pu_] && ku_.back() == ku_[nu_];
  const bool cv = kv_.front() == kv_[pv_] && kv_.back() == kv_[nv_];
  NurbsSurface s = cu ? *this : from_columns_u(*this, pu_, detail::clamp(pu_, columns_u(*this)));
  if (cv) return s;
  const NurbsSurface t = s.transposed();
  return from_columns_u(t, pv_, detail::clamp(pv_, columns_u(t))).transposed();
}

std::pair<NurbsSurface, NurbsSurface> NurbsSurface::split_u(double u) const {
  const auto [l, r] = detail::split(pu_, columns_u(*this), u);
  return {from_columns_u(*this, pu_, l), from_columns_u(*this, pu_, r)};
}

std::pair<NurbsSurface, NurbsSurface> NurbsSurface::split_v(double v) const {
  const auto [l, r] = transposed().split_u(v);
  return {l.transposed(), r.transposed()};
}

NurbsSurface NurbsSurface::sub_surface(double u0, double u1, double v0, double v1) const {
  const auto du = domain_u(), dv = domain_v();
  const auto bad = [](double a0, double a1, const std::array<double, 2>& d) {
    return !(a0 < a1) || a0 < d[0] - tol::kParametric || a1 > d[1] + tol::kParametric;
  };
  if (bad(u0, u1, du) || bad(v0, v1, dv))
    throw GeometryError("sub_surface: parameter ranges must be increasing and inside the domain");
  NurbsSurface s = clamped();
  if (u1 < du[1] - tol::kParametric) s = s.split_u(u1).first;
  if (u0 > du[0] + tol::kParametric) s = s.split_u(u0).second;
  if (v1 < dv[1] - tol::kParametric) s = s.split_v(v1).first;
  if (v0 > dv[0] + tol::kParametric) s = s.split_v(v0).second;
  return s;
}

std::pair<NurbsSurface, std::size_t> NurbsSurface::try_remove_knot_u(double u, std::size_t times,
                                                                      double tolerance) const {
  auto net = columns_u(*this);
  const std::size_t removed = detail::remove_knot(pu_, net, u, times, tolerance);
  if (removed == 0) return {*this, 0};
  return {from_columns_u(*this, pu_, net), removed};
}

std::pair<NurbsSurface, std::size_t> NurbsSurface::try_remove_knot_v(double v, std::size_t times,
                                                                      double tolerance) const {
  auto [s, removed] = transposed().try_remove_knot_u(v, times, tolerance);
  return {s.transposed(), removed};
}

NurbsSurface NurbsSurface::remove_knot_u(double u, std::size_t times, double tolerance) const {
  auto [s, removed] = try_remove_knot_u(u, times, tolerance);
  if (removed < times)
    throw GeometryError("remove_knot_u: only " + std::to_string(removed) + " of " +
                        std::to_string(times) + " removals possible within tolerance");
  return s;
}

NurbsSurface NurbsSurface::remove_knot_v(double v, std::size_t times, double tolerance) const {
  return transposed().remove_knot_u(v, times, tolerance).transposed();
}

NurbsSurface NurbsSurface::reduce_degree_u(double tolerance) const {
  return from_columns_u(*this, pu_ - 1, detail::reduce_degree(pu_, columns_u(*this), tolerance));
}

NurbsSurface NurbsSurface::reduce_degree_v(double tolerance) const {
  return transposed().reduce_degree_u(tolerance).transposed();
}

NurbsSurface NurbsSurface::refine_knots_v(const std::vector<double>& X) const {
  if (X.empty()) return *this;
  return transposed().refine_knots_u(X).transposed();
}

std::vector<NurbsSurface> NurbsSurface::decompose_bezier() const {
  const NurbsSurface s = refine_knots_u(detail::bezier_refinement_knots(pu_, ku_, nu_))
                             .refine_knots_v(detail::bezier_refinement_knots(pv_, kv_, nv_));
  const auto clamped = [](std::size_t p, double a, double b) {
    KnotVector U(2 * (p + 1), a);
    std::fill(U.begin() + static_cast<std::ptrdiff_t>(p + 1), U.end(), b);
    return U;
  };
  std::vector<NurbsSurface> patches;
  for (std::size_t ku : detail::nonempty_spans(pu_, s.ku_, s.nu_)) {
    for (std::size_t kv : detail::nonempty_spans(pv_, s.kv_, s.nv_)) {
      std::vector<Vec3> c;
      std::vector<double> w;
      for (std::size_t i = ku - pu_; i <= ku; ++i)
        for (std::size_t j = kv - pv_; j <= kv; ++j) {
          c.push_back(s.control_point(i, j));
          w.push_back(s.weight(i, j));
        }
      patches.emplace_back(pu_, pv_, clamped(pu_, s.ku_[ku], s.ku_[ku + 1]),
                           clamped(pv_, s.kv_[kv], s.kv_[kv + 1]), pu_ + 1, pv_ + 1,
                           std::move(c), std::move(w));
    }
  }
  return patches;
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
