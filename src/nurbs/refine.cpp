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

namespace {

using idx = std::ptrdiff_t;
constexpr std::size_t z(idx i) { return static_cast<std::size_t>(i); }

double binomial(idx n, idx k) {
  double r = 1.0;
  for (idx i = 1; i <= k; ++i) r = r * static_cast<double>(n - k + i) / static_cast<double>(i);
  return r;
}

// Makes knots that agree within tolerance exactly equal (keeps the first of each cluster).
KnotVector snap_knots(KnotVector U) {
  for (std::size_t i = 1; i < U.size(); ++i)
    if (U[i] - U[i - 1] < tol::kParametric) U[i] = U[i - 1];
  return U;
}

}  // namespace

RefinedPolygon clamp(std::size_t p, const KnotVector& U_in, const std::vector<Vec4>& Pw) {
  const KnotVector U = snap_knots(U_in);
  const std::size_t n_ctrl = Pw.size();
  const double a = U[p], b = U[n_ctrl];
  if (U.front() == a && U.back() == b) return {U, Pw};

  // raise the multiplicity of both domain ends to p, then cut off the outer parts
  std::vector<double> X;
  for (std::size_t c = knot_multiplicity(U, a); c < p; ++c) X.push_back(a);
  for (std::size_t c = knot_multiplicity(U, b); c < p; ++c) X.push_back(b);
  const auto r = refine_knot_vector(p, U, Pw, X);

  std::size_t e = 0;  // last index with knot == a
  while (r.knots[e] != a) ++e;
  while (r.knots[e + 1] == a) ++e;
  std::size_t f = r.knots.size() - 1;  // first index with knot == b
  while (r.knots[f] != b) --f;
  while (r.knots[f - 1] == b) --f;

  KnotVector Uc(p + 1, a);
  Uc.insert(Uc.end(), r.knots.begin() + static_cast<idx>(e + 1), r.knots.begin() + static_cast<idx>(f));
  Uc.insert(Uc.end(), p + 1, b);
  std::vector<Vec4> Qc(r.points.begin() + static_cast<idx>(e - p),
                       r.points.begin() + static_cast<idx>(f));
  return {std::move(Uc), std::move(Qc)};
}

RefinedPolygon elevate_degree(std::size_t p_in, const KnotVector& U_in,
                              const std::vector<Vec4>& Pw_in, std::size_t t_in) {
  if (t_in == 0) return {U_in, Pw_in};
  const auto [U, Pw] = clamp(p_in, U_in, Pw_in);

  // A5.9 DegreeElevateCurve
  const idx p = static_cast<idx>(p_in), t = static_cast<idx>(t_in);
  const idx n = static_cast<idx>(Pw.size()) - 1;
  const idx m = n + p + 1;
  const idx ph = p + t;
  const idx ph2 = ph / 2;

  // Bezier degree elevation coefficients
  std::vector<std::vector<double>> bezalfs(z(ph + 1), std::vector<double>(z(p + 1), 0.0));
  bezalfs[0][0] = 1.0;
  bezalfs[z(ph)][z(p)] = 1.0;
  for (idx i = 1; i <= ph2; ++i) {
    const double inv = 1.0 / binomial(ph, i);
    for (idx j = std::max<idx>(0, i - t); j <= std::min(p, i); ++j)
      bezalfs[z(i)][z(j)] = inv * binomial(p, j) * binomial(t, i - j);
  }
  for (idx i = ph2 + 1; i <= ph - 1; ++i)
    for (idx j = std::max<idx>(0, i - t); j <= std::min(p, i); ++j)
      bezalfs[z(i)][z(j)] = bezalfs[z(ph - i)][z(p - j)];

  // every Bezier segment gains t control points: (n+1)*(t+1) is an upper bound
  std::vector<Vec4> Qw(Pw.size() * z(t + 1));
  KnotVector Uh(Qw.size() + z(ph + 1));
  std::vector<Vec4> bpts(z(p + 1)), ebpts(z(ph + 1)), Nextbpts(z(std::max<idx>(p - 1, 1)));
  std::vector<double> alfs(z(std::max<idx>(p - 1, 1)));

  idx mh = ph, kind = ph + 1;
  idx r = -1, a = p, b = p + 1, cind = 1;
  double ua = U[0];
  Qw[0] = Pw[0];
  for (idx i = 0; i <= ph; ++i) Uh[z(i)] = ua;
  for (idx i = 0; i <= p; ++i) bpts[z(i)] = Pw[z(i)];

  while (b < m) {
    idx i = b;
    while (b < m && U[z(b)] == U[z(b + 1)]) ++b;
    const idx mul = b - i + 1;
    mh = mh + mul + t;
    const double ub = U[z(b)];
    const idx oldr = r;
    r = p - mul;
    const idx lbz = oldr > 0 ? (oldr + 2) / 2 : 1;
    const idx rbz = r > 0 ? ph - (r + 1) / 2 : ph;

    if (r > 0) {  // insert knot ub r times to get the Bezier segment
      const double numer = ub - ua;
      for (idx k = p; k > mul; --k) alfs[z(k - mul - 1)] = numer / (U[z(a + k)] - ua);
      for (idx j = 1; j <= r; ++j) {
        const idx save = r - j, s = mul + j;
        for (idx k = p; k >= s; --k)
          bpts[z(k)] = alfs[z(k - s)] * bpts[z(k)] + (1.0 - alfs[z(k - s)]) * bpts[z(k - 1)];
        Nextbpts[z(save)] = bpts[z(p)];
      }
    }
    for (i = lbz; i <= ph; ++i) {  // degree-elevate the Bezier segment
      ebpts[z(i)] = Vec4{};
      for (idx j = std::max<idx>(0, i - t); j <= std::min(p, i); ++j)
        ebpts[z(i)] += bezalfs[z(i)][z(j)] * bpts[z(j)];
    }
    if (oldr > 1) {  // remove knot ua oldr times
      idx first = kind - 2, last = kind;
      const double den = ub - ua;
      const double bet = (ub - Uh[z(kind - 1)]) / den;
      for (idx tr = 1; tr < oldr; ++tr) {
        i = first;
        idx j = last, kj = j - kind + 1;
        while (j - i > tr) {
          if (i < cind) {
            const double alf = (ub - Uh[z(i)]) / (ua - Uh[z(i)]);
            Qw[z(i)] = alf * Qw[z(i)] + (1.0 - alf) * Qw[z(i - 1)];
          }
          if (j >= lbz) {
            if (j - tr <= kind - ph + oldr) {
              const double gam = (ub - Uh[z(j - tr)]) / den;
              ebpts[z(kj)] = gam * ebpts[z(kj)] + (1.0 - gam) * ebpts[z(kj + 1)];
            } else {
              ebpts[z(kj)] = bet * ebpts[z(kj)] + (1.0 - bet) * ebpts[z(kj + 1)];
            }
          }
          ++i;
          --j;
          --kj;
        }
        --first;
        ++last;
      }
    }
    if (a != p)  // load the knot ua
      for (i = 0; i < ph - oldr; ++i) Uh[z(kind++)] = ua;
    for (idx j = lbz; j <= rbz; ++j) Qw[z(cind++)] = ebpts[z(j)];

    if (b < m) {  // set up for the next pass
      for (idx j = 0; j < r; ++j) bpts[z(j)] = Nextbpts[z(j)];
      for (idx j = r; j <= p; ++j) bpts[z(j)] = Pw[z(b - p + j)];
      a = b;
      ++b;
      ua = ub;
    } else {
      for (i = 0; i <= ph; ++i) Uh[z(kind + i)] = ub;
    }
  }
  const idx nh = mh - ph - 1;
  Qw.resize(z(nh + 1));
  Uh.resize(z(nh + ph + 2));
  return {std::move(Uh), std::move(Qw)};
}

Columns clamp(std::size_t p, const Columns& net) {
  Columns out;
  for (const auto& c : net.cols) {
    auto r = clamp(p, net.knots, c);
    out.knots = std::move(r.knots);
    out.cols.push_back(std::move(r.points));
  }
  return out;
}

std::pair<Columns, Columns> split(std::size_t p, const Columns& net, double u) {
  const Columns cl = clamp(p, net);
  const std::size_t n_ctrl = cl.cols.front().size();
  const double a = cl.knots[p], b = cl.knots[n_ctrl];
  if (!(u > a + tol::kParametric && u < b - tol::kParametric))
    throw GeometryError("split: parameter " + std::to_string(u) + " is not inside (" +
                        std::to_string(a) + ", " + std::to_string(b) + ")");
  for (double k : cl.knots)
    if (std::abs(k - u) < tol::kParametric) { u = k; break; }
  const std::vector<double> X(p - knot_multiplicity(cl.knots, u), u);

  Columns left, right;
  for (const auto& col : cl.cols) {
    auto r = refine_knot_vector(p, cl.knots, col, X);
    std::size_t f = 0;  // first index of u (multiplicity p)
    while (r.knots[f] != u) ++f;
    if (left.knots.empty()) {
      left.knots.assign(r.knots.begin(), r.knots.begin() + static_cast<idx>(f + p));
      left.knots.push_back(u);
      right.knots.assign(1, u);
      right.knots.insert(right.knots.end(), r.knots.begin() + static_cast<idx>(f), r.knots.end());
    }
    left.cols.emplace_back(r.points.begin(), r.points.begin() + static_cast<idx>(f));
    right.cols.emplace_back(r.points.begin() + static_cast<idx>(f - 1), r.points.end());
  }
  return {std::move(left), std::move(right)};
}

std::vector<std::size_t> nonempty_spans(std::size_t p, const KnotVector& U, std::size_t n_ctrl) {
  std::vector<std::size_t> spans;
  for (std::size_t k = p; k < n_ctrl; ++k)
    if (U[k + 1] - U[k] > tol::kParametric) spans.push_back(k);
  return spans;
}

}  // namespace nf::detail
