// Blend and match surfaces: smooth (G0/G1/G2) transitions between surface sides.
//
// Both constructions work with the "rows" of a side: the edge curve d0 and the outward
// cross-boundary derivatives d1, d2 as curves along the edge. With the side moved to v = v_min
// (and v pointing inwards), d1 = -S_v and d2 = S_vv at v_min. A surface B continuing S with
// B(u, 0) = d0, B_v = k d1, B_vv = k^2 d2 is a reparametrisation of S across the edge up to
// second order (v -> v_min - k v), hence G1 / G2.
//
// When the weights of the first rows of the control net are equal across the side, these
// curves are exact NURBS with the edge's weights and knots (the denominator does not depend on
// v up to the order used); otherwise they are interpolated from dense samples.

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

#include "nurbforge/core/error.hpp"
#include "nurbforge/core/tolerance.hpp"
#include "nurbforge/fitting/fit.hpp"
#include "nurbforge/modeling/blend.hpp"
#include "nurbforge/modeling/construct.hpp"
#include "nurbforge/nurbs/knot_vector.hpp"

namespace nf {

namespace {

bool is_u_side(SurfaceSide s) { return s == SurfaceSide::UMin || s == SurfaceSide::UMax; }
bool is_max_side(SurfaceSide s) { return s == SurfaceSide::UMax || s == SurfaceSide::VMax; }

NurbsSurface reversed_v(const NurbsSurface& s) {
  const auto& kv = s.knots_v();
  const std::size_t m = kv.size() - 1;
  KnotVector k(kv.size());
  for (std::size_t i = 0; i <= m; ++i) k[i] = kv[0] + kv[m] - kv[m - i];
  const std::size_t nu = s.num_u(), nv = s.num_v();
  std::vector<Vec3> P(nu * nv);
  std::vector<double> W(nu * nv);
  for (std::size_t i = 0; i < nu; ++i)
    for (std::size_t j = 0; j < nv; ++j) {
      P[i * nv + j] = s.control_point(i, nv - 1 - j);
      W[i * nv + j] = s.weight(i, nv - 1 - j);
    }
  return {s.degree_u(), s.degree_v(), s.knots_u(),  std::move(k),
          nu,           nv,           std::move(P), std::move(W)};
}

// The surface with the side at v = v_min and v pointing inwards; u runs along the edge.
NurbsSurface canonical(const NurbsSurface& s, SurfaceSide side) {
  NurbsSurface c = s.clamped();
  if (is_u_side(side)) c = c.transposed();
  if (is_max_side(side)) c = reversed_v(c);
  return c;
}

NurbsSurface uncanonical(const NurbsSurface& c, SurfaceSide side) {
  NurbsSurface s = is_max_side(side) ? reversed_v(c) : c;
  return is_u_side(side) ? s.transposed() : s;
}

NurbsCurve row(const NurbsSurface& c, std::size_t j) {
  std::vector<Vec3> P;
  std::vector<double> W;
  for (std::size_t i = 0; i < c.num_u(); ++i) {
    P.push_back(c.control_point(i, j));
    W.push_back(c.weight(i, j));
  }
  return {c.degree_u(), c.knots_u(), std::move(P), std::move(W)};
}

// sum_k f[k] * curve[k] for curves with identical knots and weights
NurbsCurve combine(const std::vector<const NurbsCurve*>& c, const std::vector<double>& f) {
  std::vector<Vec3> P(c.front()->num_control_points(), Vec3{0, 0, 0});
  for (std::size_t k = 0; k < c.size(); ++k)
    for (std::size_t i = 0; i < P.size(); ++i) P[i] = P[i] + f[k] * c[k]->control_points()[i];
  return {c.front()->degree(), c.front()->knots(), std::move(P), c.front()->weights()};
}

// Derivative coefficients of a clamped degree-q B-spline at its start (v-knots t):
// C' = a1 (P1 - P0), C'' = a2 ((P2 - P1) b - (P1 - P0) a1 / q) ... written out below.
struct StartDerivative {
  double c1{0.0};   // Q0 = c1 (P1 - P0)
  double c2a{0.0};  // Q1 = c2a (P2 - P1)
  double c2{0.0};   // R0 = c2 (Q1 - Q0)
};

StartDerivative start_derivative(const NurbsSurface& c) {
  const std::size_t q = c.degree_v();
  const auto& t = c.knots_v();
  StartDerivative d;
  d.c1 = static_cast<double>(q) / (t[q + 1] - t[1]);
  if (q >= 2 && c.num_v() >= 3) {
    d.c2a = static_cast<double>(q) / (t[q + 2] - t[2]);
    d.c2 = static_cast<double>(q - 1) / (t[q + 1] - t[2]);
  }
  return d;
}

// Rows {d0, d1 (outward), d2} up to order n of the canonical surface c.
std::vector<NurbsCurve> exact_rows(const NurbsSurface& c, std::size_t n) {
  std::vector<NurbsCurve> r{row(c, 0)};
  if (n == 0) return r;
  const StartDerivative d = start_derivative(c);
  const NurbsCurve p1 = row(c, 1);
  // d1 = -Q0 = -c1 (P1 - P0)
  r.push_back(combine({&r[0], &p1}, {d.c1, -d.c1}));
  if (n == 1) return r;
  if (d.c2 == 0.0) {
    r.push_back(combine({&r[0]}, {0.0}));
    return r;
  }
  const NurbsCurve p2 = row(c, 2);
  // R0 = c2 (Q1 - Q0) = c2 (c2a (P2 - P1) - c1 (P1 - P0))
  r.push_back(combine({&r[0], &p1, &p2}, {d.c2 * d.c1, -d.c2 * (d.c2a + d.c1), d.c2 * d.c2a}));
  return r;
}

bool weights_constant_across(const NurbsSurface& c, std::size_t n) {
  const std::size_t rows = std::min(n + 1, c.num_v());
  for (std::size_t i = 0; i < c.num_u(); ++i)
    for (std::size_t j = 1; j < rows; ++j)
      if (std::abs(c.weight(i, j) - c.weight(i, 0)) > 1e-12 * c.weight(i, 0)) return false;
  return true;
}

double polyline_length(const NurbsCurve& c, std::size_t n) {
  double L = 0.0;
  Vec3 prev = c.evaluate(c.domain_start());
  for (std::size_t k = 1; k <= n; ++k) {
    const double u = c.domain_start() + (c.domain_end() - c.domain_start()) *
                                            static_cast<double>(k) / static_cast<double>(n);
    const Vec3 p = c.evaluate(u);
    L += distance(prev, p);
    prev = p;
  }
  return L;
}

// Joins clamped cubic pieces that meet end to start into one curve (C0 at the joints).
NurbsCurve concatenate(const std::vector<NurbsCurve>& pieces) {
  KnotVector K(4, pieces.front().domain_start());
  std::vector<Vec3> P;
  for (std::size_t k = 0; k < pieces.size(); ++k) {
    const auto& c = pieces[k];
    const auto& kn = c.knots();
    for (std::size_t i = 4; i + 4 < kn.size(); ++i) K.push_back(kn[i]);
    for (int r = 0; r < (k + 1 == pieces.size() ? 4 : 3); ++r) K.push_back(c.domain_end());
    const auto& cp = c.control_points();
    P.insert(P.end(), cp.begin() + (k == 0 ? 0 : 1), cp.end());
  }
  return {3, std::move(K), std::move(P), {}};
}

// Rows by cubic interpolation of samples at their parameters, piecewise between the knots
// where the surface is less than C2 along the edge, refined until the midpoints are within
// tolerance (derivative rows: relative to their size over the edge length).
std::vector<NurbsCurve> sampled_rows(const NurbsSurface& c, std::size_t n, double tolerance) {
  const double v0 = c.domain_v()[0];
  const auto [kv, km] = unique_knots(c.knots_u());
  const auto p = static_cast<int>(c.degree_u());
  std::vector<std::size_t> breaks{0};
  for (std::size_t s = 1; s + 1 < kv.size(); ++s)
    if (p - km[s] < 2) breaks.push_back(s);
  breaks.push_back(kv.size() - 1);
  const auto value = [&](double u, std::size_t order) {
    const auto d = c.derivatives(u, v0, n);
    return order == 0 ? d[0][0] : order == 1 ? -1.0 * d[0][1] : d[0][2];
  };
  const double L = std::max(polyline_length(c.iso_curve_v(v0), 64), 1e-300);
  std::vector<NurbsCurve> best;
  for (std::size_t m = 4; m <= 512; m *= 2) {
    std::vector<std::vector<double>> params;  // per piece
    for (std::size_t b = 0; b + 1 < breaks.size(); ++b) {
      params.emplace_back();
      for (std::size_t s = breaks[b]; s < breaks[b + 1]; ++s)
        for (std::size_t k = 0; k < m; ++k)
          params.back().push_back(kv[s] + (kv[s + 1] - kv[s]) * static_cast<double>(k) /
                                              static_cast<double>(m));
      params.back().push_back(kv[breaks[b + 1]]);
    }
    std::vector<NurbsCurve> rows;
    bool ok = true;
    for (std::size_t o = 0; o <= n; ++o) {
      std::vector<NurbsCurve> pieces;
      double size = 0.0;
      for (const auto& pp : params) {
        std::vector<Vec3> pts;
        for (double u : pp) {
          pts.push_back(value(u, o));
          size = std::max(size, norm(pts.back()));
        }
        pieces.push_back(interpolate_curve(pts, 3, pp));
      }
      rows.push_back(pieces.size() == 1 ? pieces.front() : concatenate(pieces));
      const double tol = o == 0 ? tolerance : tolerance * std::max(size, 1e-300) / L;
      for (const auto& pp : params)
        for (std::size_t k = 0; k + 1 < pp.size() && ok; ++k) {
          const double u = 0.5 * (pp[k] + pp[k + 1]);
          if (distance(rows.back().evaluate(u), value(u, o)) > tol) ok = false;
        }
    }
    best = std::move(rows);
    if (ok) return best;
  }
  return best;
}

// Outward rows of a side, as curves along the edge in the canonical parameter.
std::vector<NurbsCurve> side_rows(const NurbsSurface& s, SurfaceSide side, std::size_t n,
                                  double tolerance) {
  const NurbsSurface c = canonical(s, side);
  return weights_constant_across(c, n) ? exact_rows(c, n) : sampled_rows(c, n, tolerance);
}

constexpr std::size_t kSamples = 64;
constexpr double kSamplesD = 64.0;

double param_at(const NurbsCurve& c, std::size_t k) {
  return c.domain_start() +
         (c.domain_end() - c.domain_start()) * static_cast<double>(k) / kSamplesD;
}

// True if b runs against a (pairing a(t) with b(1 - t) is closer).
bool opposite(const NurbsCurve& a, const NurbsCurve& b) {
  double same = 0.0, opp = 0.0;
  for (std::size_t k = 0; k <= kSamples; ++k) {
    const Vec3 p = a.evaluate(param_at(a, k));
    same += distance(p, b.evaluate(param_at(b, k)));
    opp += distance(p, b.evaluate(param_at(b, kSamples - k)));
  }
  return opp < same;
}

double mean_norm(const NurbsCurve& c) {
  double m = 0.0;
  for (std::size_t k = 0; k <= kSamples; ++k) m += norm(c.evaluate(param_at(c, k)));
  return m / (kSamplesD + 1.0);
}

void orient_and_map(std::vector<NurbsCurve>& rows, const NurbsCurve& reference, double a,
                    double b) {
  const bool flip = opposite(reference, rows[0]);
  for (auto& r : rows) r = reparametrize(flip ? r.reversed() : r, a, b);
}

void check_continuity(int continuity, double tolerance, const char* who) {
  if (continuity < 0 || continuity > 2)
    throw GeometryError(std::string(who) + ": continuity must be 0 (G0), 1 (G1) or 2 (G2)");
  if (!(tolerance > 0.0)) throw GeometryError(std::string(who) + ": tolerance must be positive");
}

// Rows of a degree-d Bezier across, from one side: P0 = d0, P1 = P0 + k/d d1,
// P2 = 2 P1 - P0 + k^2 / (d (d - 1)) d2.
std::vector<NurbsCurve> bezier_rows(const std::vector<NurbsCurve>& r, double k, std::size_t d) {
  std::vector<NurbsCurve> out{r[0]};
  const double dd = static_cast<double>(d);
  if (r.size() > 1) out.push_back(combine({&r[0], &r[1]}, {1.0, k / dd}));
  if (r.size() > 2)
    out.push_back(combine({&r[0], &r[1], &r[2]}, {1.0, 2.0 * k / dd, k * k / (dd * (dd - 1))}));
  return out;
}

}  // namespace

NurbsSurface blend_surfaces(const NurbsSurface& s0, SurfaceSide side0, const NurbsSurface& s1,
                            SurfaceSide side1, int continuity, double tension, double tolerance) {
  constexpr const char* who = "blend_surfaces";
  check_continuity(continuity, tolerance, who);
  if (!(tension > 0.0)) throw GeometryError("blend_surfaces: tension must be positive");
  const auto n = static_cast<std::size_t>(continuity);
  auto r0 = side_rows(s0, side0, n, tolerance);
  auto r1 = side_rows(s1, side1, n, tolerance);
  for (auto& r : r0) r = reparametrize(r, 0.0, 1.0);
  orient_and_map(r1, r0[0], 0.0, 1.0);

  double chord = 0.0;
  for (std::size_t k = 0; k <= kSamples; ++k) {
    const double t = static_cast<double>(k) / kSamplesD;
    chord += distance(r0[0].evaluate(t), r1[0].evaluate(t));
  }
  chord /= kSamplesD + 1.0;
  if (!(chord > 0.0)) throw GeometryError("blend_surfaces: the edges coincide");

  const std::size_t d = 2 * n + 1;
  std::vector<NurbsCurve> rows;
  std::vector<NurbsCurve> far;
  if (n == 0) {
    rows = {r0[0]};
    far = {r1[0]};
  } else {
    const double m0 = mean_norm(r0[1]), m1 = mean_norm(r1[1]);
    if (!(m0 > 0.0) || !(m1 > 0.0))
      throw GeometryError("blend_surfaces: vanishing cross-boundary derivative");
    rows = bezier_rows(r0, tension * chord / m0, d);
    far = bezier_rows(r1, tension * chord / m1, d);
  }
  rows.insert(rows.end(), far.rbegin(), far.rend());
  make_compatible(rows);

  const std::size_t nu = rows[0].num_control_points(), nv = rows.size();
  std::vector<Vec3> P(nu * nv);
  std::vector<double> W(nu * nv);
  for (std::size_t i = 0; i < nu; ++i)
    for (std::size_t j = 0; j < nv; ++j) {
      P[i * nv + j] = rows[j].control_points()[i];
      W[i * nv + j] = rows[j].weights()[i];
    }
  return {rows[0].degree(), d,           rows[0].knots(), clamped_uniform_knots(d, nv), nu, nv,
          std::move(P),     std::move(W)};
}

NurbsSurface match_surface(const NurbsSurface& surface, SurfaceSide side,
                           const NurbsSurface& target, SurfaceSide target_side, int continuity,
                           double tolerance) {
  constexpr const char* who = "match_surface";
  check_continuity(continuity, tolerance, who);
  const auto n = static_cast<std::size_t>(continuity);
  NurbsSurface c = canonical(surface, side);

  // across: degree >= n and enough rows that the opposite side keeps its derivatives
  if (c.degree_v() < std::max<std::size_t>(n, 1)) c = c.elevate_degree_v(n - c.degree_v());
  while (c.num_v() < 2 * (n + 1)) {
    const auto kv = unique_knots(c.knots_v()).first;
    std::size_t best = 0;
    for (std::size_t s = 1; s + 1 < kv.size(); ++s)
      if (kv[s + 1] - kv[s] > kv[best + 1] - kv[best]) best = s;
    c = c.refine_knots_v({0.5 * (kv[best] + kv[best + 1])});
  }
  const double v0 = c.domain_v()[0];
  const auto [a, b] = c.domain_u();

  // mean cross speed of the original surface
  double speed = 0.0;
  for (std::size_t k = 0; k <= kSamples; ++k) {
    const double u = a + (b - a) * static_cast<double>(k) / kSamplesD;
    speed += norm(c.first_derivatives(u, v0)[2]);
  }
  speed /= kSamplesD + 1.0;

  auto t = side_rows(target, target_side, n, tolerance);
  orient_and_map(t, c.iso_curve_v(v0), a, b);

  // along the edge: common degree and knots (the surface keeps its shape)
  const std::size_t p = std::max(c.degree_u(), t[0].degree());
  if (c.degree_u() < p) c = c.elevate_degree_u(p - c.degree_u());
  if (t[0].degree() < p)
    for (auto& r : t) r = r.elevate_degree(p - r.degree());
  {
    std::vector<double> vals;
    std::vector<std::size_t> mult;
    for (const auto* kv : {&c.knots_u(), &t[0].knots()}) {
      const auto [cv, cm] = unique_knots(*kv);
      for (std::size_t i = 0; i < cv.size(); ++i) {
        std::size_t k = 0;
        while (k < vals.size() && std::abs(vals[k] - cv[i]) >= tol::kParametric) ++k;
        if (k == vals.size()) {
          vals.push_back(cv[i]);
          mult.push_back(0);
        }
        mult[k] = std::max(mult[k], static_cast<std::size_t>(cm[i]));
      }
    }
    const auto missing = [&](const KnotVector& kv) {
      std::vector<double> X;
      for (std::size_t k = 0; k < vals.size(); ++k)
        for (auto have = knot_multiplicity(kv, vals[k]); have < mult[k]; ++have)
          X.push_back(vals[k]);
      std::sort(X.begin(), X.end());
      return X;
    };
    const auto xc = missing(c.knots_u());
    if (!xc.empty()) c = c.refine_knots_u(xc);
    const auto xt = missing(t[0].knots());
    for (auto& r : t) {
      if (!xt.empty()) r = r.refine_knots(xt);
      r = NurbsCurve(r.degree(), c.knots_u(), r.control_points(), r.weights());
    }
  }

  // new rows: P0 = d0, Q0 = k d1, R0 = k^2 d2 (Q, R: first / second derivative control rows)
  std::vector<NurbsCurve> rows{t[0]};
  if (n >= 1) {
    const double m = mean_norm(t[1]);
    if (!(m > 0.0) || !(speed > 0.0))
      throw GeometryError("match_surface: vanishing cross-boundary derivative");
    const double k = speed / m;
    const StartDerivative d = start_derivative(c);
    // P1 = P0 + Q0 / c1
    rows.push_back(combine({&t[0], &t[1]}, {1.0, k / d.c1}));
    if (n == 2) {
      // Q1 = Q0 + R0 / c2, P2 = P1 + Q1 / c2a
      rows.push_back(
          combine({&t[0], &t[1], &t[2]}, {1.0, k / d.c1 + k / d.c2a, k * k / (d.c2 * d.c2a)}));
    }
  }

  std::vector<Vec3> P = c.control_points();
  std::vector<double> W = c.weights();
  const std::size_t nv = c.num_v();
  for (std::size_t j = 0; j < rows.size(); ++j)
    for (std::size_t i = 0; i < c.num_u(); ++i) {
      P[i * nv + j] = rows[j].control_points()[i];
      W[i * nv + j] = rows[j].weights()[i];
    }
  const NurbsSurface out(c.degree_u(), c.degree_v(), c.knots_u(), c.knots_v(), c.num_u(), nv,
                         std::move(P), std::move(W));
  return uncanonical(out, side);
}

}  // namespace nf
