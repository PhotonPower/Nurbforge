#include "nurbforge/modeling/gordon.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <optional>
#include <string>

#include "nurbforge/analysis/closest_point.hpp"
#include "nurbforge/core/error.hpp"
#include "nurbforge/fitting/fit.hpp"
#include "nurbforge/modeling/construct.hpp"
#include "nurbforge/modeling/loft.hpp"
#include "nurbforge/nurbs/join.hpp"

namespace nf {

namespace {

NurbsCurve unit_domain(const NurbsCurve& c) { return reparametrize(c.clamped(), 0.0, 1.0); }

// Non-rational replacement within `tol`, parametrisation preserved: cubic interpolation at the
// original parameters, piecewise per knot span (rational curves are often only C1 at knots).
NurbsCurve nonrational(const NurbsCurve& c_in, double tol) {
  const NurbsCurve c = unit_domain(c_in);
  if (!c.is_rational()) return c;
  const auto breaks = unique_knots(c.knots()).first;
  std::optional<NurbsCurve> result;
  for (std::size_t s = 0; s + 1 < breaks.size(); ++s) {
    const double a = breaks[s], b = breaks[s + 1];
    for (std::size_t N = 4;; N = 2 * N - 1) {
      std::vector<double> t(N);
      std::vector<Vec3> pts(N);
      for (std::size_t j = 0; j < N; ++j) {
        t[j] = j + 1 == N ? b : a + (b - a) * static_cast<double>(j) / static_cast<double>(N - 1);
        pts[j] = c.evaluate(t[j]);
      }
      const NurbsCurve piece = interpolate_curve(pts, 3, t);
      double err = 0.0;
      for (std::size_t j = 0; j + 1 < N; ++j)
        for (double f : {0.25, 0.5, 0.75}) {
          const double tm = t[j] + f * (t[j + 1] - t[j]);
          err = std::max(err, distance(piece.evaluate(tm), c.evaluate(tm)));
        }
      if (err <= tol) {
        result = result ? join(*result, piece, tol) : piece;
        break;
      }
      if (N > 4097) throw GeometryError("gordon: cannot approximate a rational curve within tolerance");
    }
  }
  return *result;
}

// Intersection of two curves by alternating projection, started at the closest sample.
// Returns {u on a, v on b, distance}.
struct Meet {
  double u, v, gap;
};
Meet intersect(const NurbsCurve& a, const NurbsCurve& b) {
  double v = b.domain_start(), best = std::numeric_limits<double>::max();
  for (int i = 0; i <= 32; ++i) {
    const double t = b.domain_start() + (b.domain_end() - b.domain_start()) * i / 32.0;
    const double d = closest_point(a, b.evaluate(t)).distance;
    if (d < best) { best = d; v = t; }
  }
  CurvePoint pa{};
  for (int it = 0; it < 200; ++it) {
    pa = closest_point(a, b.evaluate(v));
    const CurvePoint pb = closest_point(b, pa.point);
    const bool done = std::abs(pb.u - v) <= 1e-15 * (1.0 + std::abs(v));
    v = pb.u;
    if (done || pb.distance <= 1e-14) break;
  }
  pa = closest_point(a, b.evaluate(v));
  return {pa.u, v, pa.distance};
}

double snap01(double x, double tol) {
  if (std::abs(x) <= tol) return 0.0;
  if (std::abs(x - 1.0) <= tol) return 1.0;
  return x;
}

// Brings three surfaces to common degrees and knots and returns P1 + P2 - P3.
NurbsSurface boolean_sum(NurbsSurface a, NurbsSurface b, NurbsSurface c) {
  const std::size_t pu = std::max({a.degree_u(), b.degree_u(), c.degree_u()});
  const std::size_t pv = std::max({a.degree_v(), b.degree_v(), c.degree_v()});
  for (auto* s : {&a, &b, &c}) {
    *s = s->elevate_degree_u(pu - s->degree_u()).elevate_degree_v(pv - s->degree_v());
  }
  // union of knots (per direction)
  const auto merge = [](const std::vector<const KnotVector*>& ks) {
    std::vector<double> vals;
    std::vector<std::size_t> mults;
    for (const auto* k : ks) {
      const auto [kv, km] = unique_knots(*k);
      for (std::size_t i = 0; i < kv.size(); ++i) {
        std::size_t j = 0;
        while (j < vals.size() && std::abs(vals[j] - kv[i]) >= tol::kParametric) ++j;
        if (j == vals.size()) { vals.push_back(kv[i]); mults.push_back(0); }
        mults[j] = std::max(mults[j], static_cast<std::size_t>(km[i]));
      }
    }
    return std::pair{vals, mults};
  };
  const auto missing = [](const KnotVector& k, const std::pair<std::vector<double>, std::vector<std::size_t>>& u) {
    std::vector<double> X;
    for (std::size_t j = 0; j < u.first.size(); ++j)
      for (auto have = knot_multiplicity(k, u.first[j]); have < u.second[j]; ++have) X.push_back(u.first[j]);
    return X;
  };
  const auto uu = merge({&a.knots_u(), &b.knots_u(), &c.knots_u()});
  const auto uv = merge({&a.knots_v(), &b.knots_v(), &c.knots_v()});
  for (auto* s : {&a, &b, &c}) *s = s->refine_knots_u(missing(s->knots_u(), uu)).refine_knots_v(missing(s->knots_v(), uv));
  std::vector<Vec3> P(a.control_points().size());
  for (std::size_t i = 0; i < P.size(); ++i)
    P[i] = a.control_points()[i] + b.control_points()[i] - c.control_points()[i];
  return {pu, pv, a.knots_u(), a.knots_v(), a.num_u(), a.num_v(), std::move(P)};
}

}  // namespace

NurbsSurface gordon(const std::vector<NurbsCurve>& u_in, const std::vector<NurbsCurve>& v_in,
                    double tolerance) {
  if (u_in.size() < 2 || v_in.size() < 2)
    throw GeometryError("gordon: need at least two curves in each direction");
  std::vector<NurbsCurve> U, V;
  for (const auto& c : u_in) U.push_back(nonrational(c, tolerance));
  for (const auto& c : v_in) V.push_back(nonrational(c, tolerance));
  const std::size_t K = U.size(), L = V.size();

  // network parameters: u_l from the first u-curve, v_k from the first v-curve
  std::vector<double> up(L), vp(K);
  for (std::size_t l = 0; l < L; ++l) {
    const Meet m = intersect(U[0], V[l]);
    if (m.gap > tolerance) throw GeometryError("gordon: v-curve " + std::to_string(l) + " does not meet the first u-curve");
    up[l] = m.u;
  }
  for (std::size_t k = 0; k < K; ++k) {
    const Meet m = intersect(V[0], U[k]);
    if (m.gap > tolerance) throw GeometryError("gordon: u-curve " + std::to_string(k) + " does not meet the first v-curve");
    vp[k] = m.u;
  }
  // order the curves by their network parameters
  std::vector<std::size_t> ou(L), ov(K);
  std::iota(ou.begin(), ou.end(), 0);
  std::iota(ov.begin(), ov.end(), 0);
  std::sort(ou.begin(), ou.end(), [&](std::size_t x, std::size_t y) { return up[x] < up[y]; });
  std::sort(ov.begin(), ov.end(), [&](std::size_t x, std::size_t y) { return vp[x] < vp[y]; });
  std::vector<NurbsCurve> Us, Vs;
  std::vector<double> u_params, v_params;
  for (auto l : ou) { Vs.push_back(V[l]); u_params.push_back(snap01(up[l], 1e-9)); }
  for (auto k : ov) { Us.push_back(U[k]); v_params.push_back(snap01(vp[k], 1e-9)); }
  if (u_params.front() != 0.0 || u_params.back() != 1.0 || v_params.front() != 0.0 || v_params.back() != 1.0)
    throw GeometryError("gordon: the outermost curves must lie on the ends of the crossing curves");
  for (std::size_t i = 1; i < L; ++i)
    if (!(u_params[i] > u_params[i - 1])) throw GeometryError("gordon: two v-curves meet the u-curves at the same parameter");
  for (std::size_t i = 1; i < K; ++i)
    if (!(v_params[i] > v_params[i - 1])) throw GeometryError("gordon: two u-curves meet the v-curves at the same parameter");

  // consistency of the network: C_k(u_l) = D_l(v_k)
  std::vector<std::vector<Vec3>> Q(K, std::vector<Vec3>(L));
  for (std::size_t k = 0; k < K; ++k)
    for (std::size_t l = 0; l < L; ++l) {
      const Vec3 a = Us[k].evaluate(u_params[l]), b = Vs[l].evaluate(v_params[k]);
      if (distance(a, b) > tolerance)
        throw GeometryError("gordon: inconsistent network – u-curve " + std::to_string(k) + " at u = " +
                            std::to_string(u_params[l]) + " misses v-curve " + std::to_string(l) +
                            " at v = " + std::to_string(v_params[k]) + " (the curves must share "
                            "their intersection parameters)");
      Q[k][l] = 0.5 * (a + b);
    }

  const std::size_t p = std::min<std::size_t>(3, L - 1), q = std::min<std::size_t>(3, K - 1);
  const NurbsSurface L1 = loft(Us, q, v_params);
  const NurbsSurface L2 = loft(Vs, p, u_params).transposed();
  std::vector<NurbsCurve> rows;
  for (std::size_t k = 0; k < K; ++k) rows.push_back(interpolate_curve(Q[k], p, u_params));
  const NurbsSurface T = loft(rows, q, v_params);
  return boolean_sum(L1, L2, T);
}

NurbsSurface coons(const NurbsCurve& c0, const NurbsCurve& c1_in, const NurbsCurve& d0_in,
                   const NurbsCurve& d1_in, double tolerance) {
  const auto start = [](const NurbsCurve& c) { return c.evaluate(c.domain_start()); };
  const auto end = [](const NurbsCurve& c) { return c.evaluate(c.domain_end()); };
  NurbsCurve c1 = c1_in, d0 = d0_in, d1 = d1_in;
  // c1 runs like c0 (decided by the opposite corners)
  if (distance(start(c0), start(c1)) + distance(end(c0), end(c1)) >
      distance(start(c0), end(c1)) + distance(end(c0), start(c1)))
    c1 = c1.reversed();
  // d0 from c0's start to c1's start, d1 from c0's end to c1's end
  if (distance(start(d0), start(c0)) + distance(end(d0), start(c1)) >
      distance(end(d0), start(c0)) + distance(start(d0), start(c1)))
    d0 = d0.reversed();
  if (distance(start(d1), end(c0)) + distance(end(d1), end(c1)) >
      distance(end(d1), end(c0)) + distance(start(d1), end(c1)))
    d1 = d1.reversed();
  const double gap = std::max({distance(start(d0), start(c0)), distance(end(d0), start(c1)),
                               distance(start(d1), end(c0)), distance(end(d1), end(c1))});
  if (gap > tolerance)
    throw GeometryError("coons: the boundary curves do not form a closed loop (corner gap " +
                        std::to_string(gap) + ")");
  return gordon({c0, c1}, {d0, d1}, tolerance);
}

}  // namespace nf
