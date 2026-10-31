#include "nurbforge/nurbs/curve.hpp"

#include <algorithm>
#include <cmath>

#include "nurbforge/core/error.hpp"
#include "nurbforge/geometry/transform.hpp"
#include "nurbforge/nurbs/basis.hpp"
#include "refine_detail.hpp"

namespace nf {

namespace {
double binomial(std::size_t n, std::size_t k) {
  double r = 1.0;
  for (std::size_t i = 1; i <= k; ++i) r = r * static_cast<double>(n - k + i) / static_cast<double>(i);
  return r;
}
}  // namespace

NurbsCurve::NurbsCurve(std::size_t degree, KnotVector knots, std::vector<Vec3> control_points,
                       std::vector<double> weights)
    : degree_(degree), knots_(std::move(knots)), ctrl_(std::move(control_points)),
      weights_(std::move(weights)) {
  if (weights_.empty()) weights_.assign(ctrl_.size(), 1.0);
  if (weights_.size() != ctrl_.size())
    throw GeometryError("number of weights must match number of control points");
  for (double w : weights_)
    if (!(w > 0.0)) throw GeometryError("weights must be strictly positive");
  validate_knots(knots_, degree_, ctrl_.size());
}

NurbsCurve NurbsCurve::from_points(std::size_t degree, std::vector<Vec3> control_points) {
  auto knots = clamped_uniform_knots(degree, control_points.size());
  return {degree, std::move(knots), std::move(control_points)};
}

bool NurbsCurve::is_rational() const noexcept {
  // Weights are scale-invariant (w and c*w give the same curve), so compare relatively.
  const double w0 = weights_.front();
  return std::any_of(weights_.begin(), weights_.end(),
                     [&](double w) { return std::abs(w - w0) > tol::kLength * std::abs(w0); });
}

std::vector<Vec4> NurbsCurve::homogeneous_control_points() const {
  std::vector<Vec4> pw(ctrl_.size());
  for (std::size_t i = 0; i < ctrl_.size(); ++i) pw[i] = Vec4::from_point(ctrl_[i], weights_[i]);
  return pw;
}

NurbsCurve NurbsCurve::from_homogeneous(std::size_t degree, KnotVector knots,
                                        const std::vector<Vec4>& pw) {
  std::vector<Vec3> pts(pw.size());
  std::vector<double> w(pw.size());
  for (std::size_t i = 0; i < pw.size(); ++i) { pts[i] = pw[i].project(); w[i] = pw[i].w; }
  return {degree, std::move(knots), std::move(pts), std::move(w)};
}

Vec3 NurbsCurve::evaluate(double u) const {
  const std::size_t n = ctrl_.size() - 1;
  const std::size_t span = find_span(n, degree_, u, knots_);
  const auto N = basis_functions(span, u, degree_, knots_);
  Vec4 cw;
  for (std::size_t j = 0; j <= degree_; ++j) {
    const std::size_t i = span - degree_ + j;
    cw += N[j] * Vec4::from_point(ctrl_[i], weights_[i]);
  }
  return cw.project();
}

std::vector<Vec3> NurbsCurve::derivatives(double u, std::size_t order) const {
  const std::size_t n = ctrl_.size() - 1;
  const std::size_t span = find_span(n, degree_, u, knots_);
  const auto nders = basis_function_derivatives(span, u, degree_, order, knots_);

  // derivatives of the homogeneous curve: A^(k) (xyz part) and w^(k)
  std::vector<Vec3> A(order + 1);
  std::vector<double> wd(order + 1, 0.0);
  for (std::size_t k = 0; k <= order; ++k) {
    for (std::size_t j = 0; j <= degree_; ++j) {
      const std::size_t i = span - degree_ + j;
      A[k] += nders[k][j] * weights_[i] * ctrl_[i];
      wd[k] += nders[k][j] * weights_[i];
    }
  }
  // A4.2: rational derivatives
  std::vector<Vec3> ck(order + 1);
  for (std::size_t k = 0; k <= order; ++k) {
    Vec3 v = A[k];
    for (std::size_t i = 1; i <= k; ++i) v -= binomial(k, i) * wd[i] * ck[k - i];
    ck[k] = v / wd[0];
  }
  return ck;
}

Vec3 NurbsCurve::tangent(double u) const { return normalized(derivatives(u, 1)[1]); }

std::vector<Vec3> NurbsCurve::sample(std::size_t n) const {
  if (n < 2) throw GeometryError("sample() needs n >= 2");
  std::vector<Vec3> pts(n);
  const double a = domain_start(), b = domain_end();
  for (std::size_t i = 0; i < n; ++i)
    pts[i] = evaluate(a + (b - a) * static_cast<double>(i) / static_cast<double>(n - 1));
  return pts;
}

NurbsCurve NurbsCurve::insert_knot(double u, std::size_t r) const {
  if (u < domain_start() || u > domain_end()) throw GeometryError("knot outside curve domain");
  const std::size_t p = degree_;
  const std::size_t s = knot_multiplicity(knots_, u);
  if (r == 0) return *this;
  if (r + s > p) throw GeometryError("knot multiplicity would exceed degree");

  const std::size_t np = ctrl_.size() - 1;
  const std::size_t k = find_span(np, p, u, knots_);
  const auto Pw = homogeneous_control_points();
  const std::size_t mp = np + p + 1;

  KnotVector UQ(knots_.size() + r);
  for (std::size_t i = 0; i <= k; ++i) UQ[i] = knots_[i];
  for (std::size_t i = 1; i <= r; ++i) UQ[k + i] = u;
  for (std::size_t i = k + 1; i <= mp; ++i) UQ[i + r] = knots_[i];

  std::vector<Vec4> Qw(Pw.size() + r);
  for (std::size_t i = 0; i <= k - p; ++i) Qw[i] = Pw[i];
  for (std::size_t i = k - s; i <= np; ++i) Qw[i + r] = Pw[i];
  std::vector<Vec4> Rw(p - s + 1);
  for (std::size_t i = 0; i <= p - s; ++i) Rw[i] = Pw[k - p + i];

  std::size_t L = 0;
  for (std::size_t j = 1; j <= r; ++j) {
    L = k - p + j;
    for (std::size_t i = 0; i <= p - j - s; ++i) {
      const double alpha = (u - knots_[L + i]) / (knots_[i + k + 1] - knots_[L + i]);
      Rw[i] = alpha * Rw[i + 1] + (1.0 - alpha) * Rw[i];
    }
    Qw[L] = Rw[0];
    Qw[k + r - j - s] = Rw[p - j - s];
  }
  for (std::size_t i = L + 1; i < k - s; ++i) Qw[i] = Rw[i - L];

  return from_homogeneous(p, std::move(UQ), Qw);
}

NurbsCurve NurbsCurve::refine_knots(const std::vector<double>& X) const {
  auto r = detail::refine_knot_vector(degree_, knots_, homogeneous_control_points(), X);
  return from_homogeneous(degree_, std::move(r.knots), r.points);
}

std::vector<NurbsCurve> NurbsCurve::decompose_bezier() const {
  const std::size_t p = degree_;
  const auto X = detail::bezier_refinement_knots(p, knots_, ctrl_.size());
  const auto r = detail::refine_knot_vector(p, knots_, homogeneous_control_points(), X);
  std::vector<NurbsCurve> segments;
  for (std::size_t k : detail::nonempty_spans(p, r.knots, r.points.size())) {
    KnotVector U(2 * (p + 1), r.knots[k]);
    std::fill(U.begin() + static_cast<std::ptrdiff_t>(p + 1), U.end(), r.knots[k + 1]);
    std::vector<Vec4> pw(r.points.begin() + static_cast<std::ptrdiff_t>(k - p),
                         r.points.begin() + static_cast<std::ptrdiff_t>(k + 1));
    segments.push_back(from_homogeneous(p, std::move(U), pw));
  }
  return segments;
}

NurbsCurve NurbsCurve::reversed() const {
  const double a = domain_start(), b = domain_end();
  KnotVector k(knots_.size());
  for (std::size_t i = 0; i < knots_.size(); ++i) k[i] = a + b - knots_[knots_.size() - 1 - i];
  std::vector<Vec3> c(ctrl_.rbegin(), ctrl_.rend());
  std::vector<double> w(weights_.rbegin(), weights_.rend());
  return {degree_, std::move(k), std::move(c), std::move(w)};
}

NurbsCurve NurbsCurve::transformed(const Transform& t) const {
  std::vector<Vec3> c(ctrl_.size());
  std::transform(ctrl_.begin(), ctrl_.end(), c.begin(), [&](const Vec3& p) { return t.apply(p); });
  return {degree_, knots_, std::move(c), weights_};
}

}  // namespace nf
