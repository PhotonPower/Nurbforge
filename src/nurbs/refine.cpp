#include <algorithm>
#include <cmath>
#include <cstddef>
#include <string>

#include "nurbforge/core/error.hpp"
#include "refine_detail.hpp"

namespace nf::detail {

RefinedPolygon refine_knot_vector(std::size_t p, const KnotVector& U, const std::vector<Vec4>& Pw,
                                  std::vector<double> X) {
  if (X.empty()) return {U, Pw};
  std::sort(X.begin(), X.end());

  const std::size_t n_ctrl = Pw.size();
  const double a_dom = U[p], b_dom = U[n_ctrl];
  for (double x : X)
    if (x < a_dom - tol::kParametric || x > b_dom + tol::kParametric)
      throw GeometryError("refinement knot " + std::to_string(x) + " outside domain [" +
                          std::to_string(a_dom) + ", " + std::to_string(b_dom) + "]");
  // Snap to existing knots / domain ends so equal values are exactly equal, then check
  // multiplicities: every value may end up with at most p copies.
  for (double& x : X) {
    x = std::clamp(x, a_dom, b_dom);
    for (double k : U)
      if (std::abs(k - x) < tol::kParametric) { x = k; break; }
  }
  for (std::size_t i = 1; i < X.size(); ++i)
    if (X[i] - X[i - 1] < tol::kParametric) X[i] = X[i - 1];
  for (std::size_t i = 0; i < X.size();) {
    std::size_t j = i;
    while (j < X.size() && X[j] == X[i]) ++j;
    const std::size_t total = knot_multiplicity(U, X[i]) + (j - i);
    if (total > p)
      throw GeometryError("refinement would raise multiplicity of knot " + std::to_string(X[i]) +
                          " to " + std::to_string(total) + " (degree " + std::to_string(p) + ")");
    i = j;
  }

  // A5.4 RefineKnotVectCurve, written with signed indices.
  using idx = std::ptrdiff_t;
  const idx pp = static_cast<idx>(p);
  const idx n = static_cast<idx>(n_ctrl) - 1;
  const idx m = n + pp + 1;
  const idx r = static_cast<idx>(X.size()) - 1;
  const auto at = [](const auto& v, idx i) -> const auto& { return v[static_cast<std::size_t>(i)]; };

  const idx a = static_cast<idx>(find_span(n_ctrl - 1, p, X.front(), U));
  const idx b = static_cast<idx>(find_span(n_ctrl - 1, p, X.back(), U)) + 1;

  KnotVector Ub(U.size() + X.size());
  std::vector<Vec4> Qw(Pw.size() + X.size());
  auto ub = [&](idx i) -> double& { return Ub[static_cast<std::size_t>(i)]; };
  auto qw = [&](idx i) -> Vec4& { return Qw[static_cast<std::size_t>(i)]; };

  for (idx j = 0; j <= a - pp; ++j) qw(j) = at(Pw, j);
  for (idx j = b - 1; j <= n; ++j) qw(j + r + 1) = at(Pw, j);
  for (idx j = 0; j <= a; ++j) ub(j) = at(U, j);
  for (idx j = b + pp; j <= m; ++j) ub(j + r + 1) = at(U, j);

  idx i = b + pp - 1;
  idx k = b + pp + r;
  for (idx j = r; j >= 0; --j) {
    const double x = at(X, j);
    while (x <= at(U, i) && i > a) {
      qw(k - pp - 1) = at(Pw, i - pp - 1);
      ub(k) = at(U, i);
      --k;
      --i;
    }
    qw(k - pp - 1) = qw(k - pp);
    for (idx l = 1; l <= pp; ++l) {
      const idx ind = k - pp + l;
      double alpha = ub(k + l) - x;
      if (std::abs(alpha) == 0.0) {
        qw(ind - 1) = qw(ind);
      } else {
        alpha /= ub(k + l) - at(U, i - pp + l);
        qw(ind - 1) = alpha * qw(ind - 1) + (1.0 - alpha) * qw(ind);
      }
    }
    ub(k) = x;
    --k;
  }
  return {std::move(Ub), std::move(Qw)};
}

std::vector<double> bezier_refinement_knots(std::size_t p, const KnotVector& U,
                                            std::size_t n_ctrl) {
  const double a = U[p], b = U[n_ctrl];
  std::vector<double> X;
  const auto [values, mults] = unique_knots(U);
  for (std::size_t i = 0; i < values.size(); ++i) {
    const double v = values[i];
    if (v < a || v > b) continue;
    const auto mult = static_cast<std::size_t>(mults[i]);
    for (std::size_t c = mult; c < p; ++c) X.push_back(v);
  }
  return X;
}

std::vector<std::size_t> nonempty_spans(std::size_t p, const KnotVector& U, std::size_t n_ctrl) {
  std::vector<std::size_t> spans;
  for (std::size_t k = p; k < n_ctrl; ++k)
    if (U[k + 1] - U[k] > tol::kParametric) spans.push_back(k);
  return spans;
}

}  // namespace nf::detail
