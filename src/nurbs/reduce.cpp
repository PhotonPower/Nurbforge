// Knot removal (A5.8) and degree reduction for one or several control polygons.
#include <algorithm>
#include <cmath>
#include <limits>
#include <string>

#include "nurbforge/core/error.hpp"
#include "refine_detail.hpp"

namespace nf::detail {

namespace {

using idx = std::ptrdiff_t;
constexpr std::size_t z(idx i) { return static_cast<std::size_t>(i); }

double dist4(const Vec4& a, const Vec4& b) {
  const double dx = a.x - b.x, dy = a.y - b.y, dz = a.z - b.z, dw = a.w - b.w;
  return std::sqrt(dx * dx + dy * dy + dz * dz + dw * dw);
}

Vec4 lin(double s, const Vec4& a, double t, const Vec4& b) { return s * a + t * b; }

// Scale factor between homogeneous and Cartesian deviations: d4 = d * factor.
// Non-rational (all weights equal to w): homogeneous distances are exactly w times the Cartesian
// ones, and the control-point deviation bounds the curve deviation. Rational: P&T Eq. 5.30,
// factor = w_min / (1 + max |P_i|).
double homogeneous_factor(const std::vector<std::vector<Vec4>>& cols) {
  double wmin = std::numeric_limits<double>::max(), wmax = 0.0, pmax = 0.0;
  for (const auto& c : cols)
    for (const auto& q : c) {
      wmin = std::min(wmin, q.w);
      wmax = std::max(wmax, q.w);
      pmax = std::max(pmax, norm(q.project()));
    }
  if (wmax - wmin <= tol::kLength * wmax) return wmin;
  return wmin / (1.0 + pmax);
}

double to_homogeneous(double d, const std::vector<std::vector<Vec4>>& cols) {
  return d * homogeneous_factor(cols);
}
double to_cartesian(double d4, const std::vector<std::vector<Vec4>>& cols) {
  return d4 / homogeneous_factor(cols);
}

// A5.8 for one polygon (knots are not modified). r: last index of u, s: multiplicity.
// Returns the number of removals t and overwrites Pw with the first n+1-t new points.
std::size_t remove_knot_polygon(idx p, const KnotVector& U, std::vector<Vec4>& Pw, double u,
                                idx r, idx s, idx num, double TOL) {
  const idx n = static_cast<idx>(Pw.size()) - 1;
  const idx ord = p + 1;
  const idx fout = (2 * r - s - p) / 2;
  idx last = r - s, first = r - p;
  std::vector<Vec4> temp(z(2 * p + 1));
  idx t = 0;
  for (; t < num; ++t) {
    const idx off = first - 1;
    temp[0] = Pw[z(off)];
    temp[z(last + 1 - off)] = Pw[z(last + 1)];
    idx i = first, j = last, ii = 1, jj = last - off;
    while (j - i > t) {
      const double alfi = (u - U[z(i)]) / (U[z(i + ord + t)] - U[z(i)]);
      const double alfj = (u - U[z(j - t)]) / (U[z(j + ord)] - U[z(j - t)]);
      temp[z(ii)] = (1.0 / alfi) * lin(1.0, Pw[z(i)], -(1.0 - alfi), temp[z(ii - 1)]);
      temp[z(jj)] = (1.0 / (1.0 - alfj)) * lin(1.0, Pw[z(j)], -alfj, temp[z(jj + 1)]);
      ++i; ++ii;
      --j; --jj;
    }
    bool removable = false;
    if (j - i < t) {
      removable = dist4(temp[z(ii - 1)], temp[z(jj + 1)]) <= TOL;
    } else {
      const double alfi = (u - U[z(i)]) / (U[z(i + ord + t)] - U[z(i)]);
      removable = dist4(Pw[z(i)], lin(alfi, temp[z(ii + t + 1)], 1.0 - alfi, temp[z(ii - 1)])) <= TOL;
    }
    if (!removable) break;
    i = first;
    j = last;
    while (j - i > t) {
      Pw[z(i)] = temp[z(i - off)];
      Pw[z(j)] = temp[z(j - off)];
      ++i;
      --j;
    }
    --first;
    ++last;
  }
  if (t == 0) return 0;
  idx j = fout, i = fout;
  for (idx k = 1; k < t; ++k) {
    if (k % 2 == 1) ++i; else --j;
  }
  for (idx k = i + 1; k <= n; ++k) Pw[z(j++)] = Pw[z(k)];
  Pw.resize(Pw.size() - z(t));
  return z(t);
}

// Degree reduction of one Bezier segment Q_0..Q_p -> P_0..P_{p-1}, endpoints kept.
std::vector<Vec4> reduce_bezier(const std::vector<Vec4>& Q) {
  const idx p = static_cast<idx>(Q.size()) - 1;
  const idx q = p - 1;
  const auto alpha = [p](idx i) { return static_cast<double>(i) / static_cast<double>(p); };
  std::vector<Vec4> P(z(q + 1));
  P[0] = Q[0];
  P[z(q)] = Q[z(p)];
  const idx r = q / 2;  // pivot
  for (idx i = 1; i < r; ++i)  // from the left: Q_i = a_i P_{i-1} + (1-a_i) P_i
    P[z(i)] = (1.0 / (1.0 - alpha(i))) * lin(1.0, Q[z(i)], -alpha(i), P[z(i - 1)]);
  for (idx i = q - 1; i > r; --i)  // from the right: Q_{i+1} = a P_i + (1-a) P_{i+1}
    P[z(i)] = (1.0 / alpha(i + 1)) * lin(1.0, Q[z(i + 1)], -(1.0 - alpha(i + 1)), P[z(i + 1)]);
  if (r > 0 && r < q) {  // pivot: average of left and right estimate
    const Vec4 L = (1.0 / (1.0 - alpha(r))) * lin(1.0, Q[z(r)], -alpha(r), P[z(r - 1)]);
    const Vec4 R = (1.0 / alpha(r + 1)) * lin(1.0, Q[z(r + 1)], -(1.0 - alpha(r + 1)), P[z(r + 1)]);
    P[z(r)] = lin(0.5, L, 0.5, R);
  }
  return P;
}

// Refines both polygons (same degree, clamped, same domain) to the union of their knots and
// returns the maximum homogeneous control-point distance.
double polygon_distance(std::size_t p, const KnotVector& Ua, const std::vector<Vec4>& Pa,
                        const KnotVector& Ub, const std::vector<Vec4>& Pb) {
  const auto missing = [](const KnotVector& from, const KnotVector& into) {
    std::vector<double> X;
    const auto [vals, mults] = unique_knots(from);
    for (std::size_t i = 0; i < vals.size(); ++i) {
      const auto have = knot_multiplicity(into, vals[i]);
      for (std::size_t c = have; c < static_cast<std::size_t>(mults[i]); ++c) X.push_back(vals[i]);
    }
    return X;
  };
  const auto a = refine_knot_vector(p, Ua, Pa, missing(Ub, Ua));
  const auto b = refine_knot_vector(p, Ub, Pb, missing(Ua, Ub));
  if (a.points.size() != b.points.size()) return std::numeric_limits<double>::infinity();
  double d = 0.0;
  for (std::size_t i = 0; i < a.points.size(); ++i) d = std::max(d, dist4(a.points[i], b.points[i]));
  return d;
}

}  // namespace

std::size_t remove_knot(std::size_t p, Columns& net, double u, std::size_t times,
                        double tolerance) {
  KnotVector& U = net.knots;
  const std::size_t n_ctrl = net.cols.front().size();
  const double a = U[p], b = U[n_ctrl];
  if (!(u > a + tol::kParametric && u < b - tol::kParametric))
    throw GeometryError("remove_knot: " + std::to_string(u) + " is not an interior knot");
  idx r = -1;
  for (std::size_t k = 0; k < U.size(); ++k)
    if (std::abs(U[k] - u) < tol::kParametric) r = static_cast<idx>(k);
  if (r < 0) throw GeometryError("remove_knot: " + std::to_string(u) + " is not a knot");
  u = U[z(r)];
  const idx s = static_cast<idx>(knot_multiplicity(U, u));
  const idx num = std::min(static_cast<idx>(times), s);
  if (num == 0) return 0;

  const double TOL = to_homogeneous(tolerance, net.cols);
  const idx pp = static_cast<idx>(p);
  // the removal count must be common to all columns: dry run first
  std::size_t t = static_cast<std::size_t>(num);
  for (const auto& c : net.cols) {
    auto copy = c;
    t = std::min(t, remove_knot_polygon(pp, U, copy, u, r, s, num, TOL));
    if (t == 0) return 0;
  }
  for (auto& c : net.cols) remove_knot_polygon(pp, U, c, u, r, s, static_cast<idx>(t), TOL);
  U.erase(U.begin() + (r - static_cast<idx>(t) + 1), U.begin() + r + 1);
  return t;
}

Columns reduce_degree(std::size_t p, const Columns& net_in, double tolerance) {
  if (p < 2) throw GeometryError("reduce_degree: degree must be >= 2");
  const std::size_t q = p - 1;

  // clamped originals (common knots for all columns)
  Columns orig;
  for (const auto& c : net_in.cols) {
    auto r = clamp(p, net_in.knots, c);
    orig.knots = std::move(r.knots);
    orig.cols.push_back(std::move(r.points));
  }
  const KnotVector& U = orig.knots;

  // Bezier-wise reduction -> C0 B-spline of degree q
  const auto X = bezier_refinement_knots(p, U, orig.cols.front().size());
  Columns c0;
  for (const auto& col : orig.cols) {
    const auto rf = refine_knot_vector(p, U, col, X);
    const auto spans = nonempty_spans(p, rf.knots, rf.points.size());
    std::vector<Vec4> pts;
    if (c0.knots.empty()) {
      c0.knots.assign(q + 1, rf.knots[spans.front()]);
      for (std::size_t k = 1; k < spans.size(); ++k) c0.knots.insert(c0.knots.end(), q, rf.knots[spans[k]]);
      c0.knots.insert(c0.knots.end(), q + 1, rf.knots[spans.back() + 1]);
    }
    for (std::size_t k : spans) {
      const std::vector<Vec4> Q(rf.points.begin() + static_cast<idx>(k - p),
                                rf.points.begin() + static_cast<idx>(k + 1));
      const auto P = reduce_bezier(Q);
      pts.insert(pts.end(), P.begin() + (pts.empty() ? 0 : 1), P.end());
    }
    c0.cols.push_back(std::move(pts));
  }

  // rigorous error: elevate back and compare control polygons
  const auto within = [&](const Columns& red) {
    double err4 = 0.0;
    for (std::size_t i = 0; i < red.cols.size(); ++i) {
      const auto e = elevate_degree(q, red.knots, red.cols[i], 1);
      err4 = std::max(err4, polygon_distance(p, e.knots, e.points, U, orig.cols[i]));
    }
    return to_cartesian(err4, orig.cols) <= tolerance;
  };
  if (!within(c0))
    throw GeometryError("reduce_degree: curve cannot be degree-reduced within tolerance " +
                        std::to_string(tolerance));

  // restore the original continuity C^(p-s): multiplicity s -> s-1 in the result. Simple knots
  // (s = 1, impossible for an elevated curve) are kept with multiplicity 1, i.e. C^(p-2).
  Columns smooth = c0;
  const auto [vals, mults] = unique_knots(U);
  for (std::size_t i = 1; i + 1 < vals.size(); ++i) {
    const auto target = std::max<std::size_t>(static_cast<std::size_t>(mults[i]) - 1, 1);
    if (target < q) remove_knot(q, smooth, vals[i], q - target, tolerance);
  }
  return within(smooth) ? smooth : c0;
}

}  // namespace nf::detail
