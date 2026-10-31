#include "nurbforge/analysis/properties.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <string>
#include <vector>

#include "nurbforge/core/error.hpp"

namespace nf {

namespace {

// ------------------------------------------------------------------ Gauss-Kronrod 7/15
// Nodes and weights from QUADPACK (qk15). Kronrod nodes xgk[0..7] (xgk[7] = 0); the Gauss
// 7-point rule uses the odd-indexed nodes xgk[1], xgk[3], xgk[5] and the centre.
constexpr std::array<double, 8> xgk{
    0.991455371120812639206854697526329, 0.949107912342758524526189684047851,
    0.864864423359769072789712788640926, 0.741531185599394439863864773280788,
    0.586087235467691130294144845693013, 0.405845151377397166906606412076961,
    0.207784955007898467600689403773245, 0.0};
constexpr std::array<double, 8> wgk{
    0.022935322010529224963732008058970, 0.063092092629978553290700663189204,
    0.104790010322250183839876322541518, 0.140653259715525918745189590510238,
    0.169004726639267902826583426598550, 0.190350578064785409913256402421014,
    0.204432940075298892414161999234649, 0.209482141084727828012999174891714};
constexpr std::array<double, 4> wg{
    0.129484966168869693270611432679082, 0.279705391489276667901467771423780,
    0.381830050505118944950369775488975, 0.417959183673469387755102040816327};

// The 15 nodes on [-1, 1] with Kronrod and (embedded) Gauss weights.
struct Rule15 {
  std::array<double, 15> x{}, wk{}, wgauss{};
  constexpr Rule15() {
    for (std::size_t j = 0; j < 7; ++j) {
      x[j] = -xgk[j];
      x[14 - j] = xgk[j];
      wk[j] = wk[14 - j] = wgk[j];
      wgauss[j] = wgauss[14 - j] = (j % 2 == 1) ? wg[j / 2] : 0.0;
    }
    x[7] = 0.0;
    wk[7] = wgk[7];
    wgauss[7] = wg[3];
  }
};
constexpr Rule15 kRule{};

template <class F>
double gk15(const F& f, double a, double b, double& err) {
  const double c = 0.5 * (a + b), h = 0.5 * (b - a);
  double K = 0.0, G = 0.0;
  for (std::size_t i = 0; i < 15; ++i) {
    const double fx = f(c + h * kRule.x[i]);
    K += kRule.wk[i] * fx;
    G += kRule.wgauss[i] * fx;
  }
  err = std::abs(h * (K - G));
  return h * K;
}

// Adaptive bisection; `density` is the accepted error per unit parameter length.
template <class F>
double integrate_1d(const F& f, double a, double b, double density, int depth = 0) {
  double err = 0.0;
  const double K = gk15(f, a, b, err);
  if (err <= density * std::abs(b - a) || depth >= 50) return K;
  const double m = 0.5 * (a + b);
  return integrate_1d(f, a, m, density, depth + 1) + integrate_1d(f, m, b, density, depth + 1);
}

template <class F>
double gk15_2d(const F& f, double ua, double ub, double va, double vb, double& err) {
  const double cu = 0.5 * (ua + ub), hu = 0.5 * (ub - ua);
  const double cv = 0.5 * (va + vb), hv = 0.5 * (vb - va);
  double K = 0.0, G = 0.0;
  for (std::size_t i = 0; i < 15; ++i)
    for (std::size_t j = 0; j < 15; ++j) {
      const double fx = f(cu + hu * kRule.x[i], cv + hv * kRule.x[j]);
      K += kRule.wk[i] * kRule.wk[j] * fx;
      G += kRule.wgauss[i] * kRule.wgauss[j] * fx;
    }
  err = std::abs(hu * hv * (K - G));
  return hu * hv * K;
}

template <class F>
double integrate_2d(const F& f, double ua, double ub, double va, double vb, double density,
                    int depth = 0) {
  double err = 0.0;
  const double K = gk15_2d(f, ua, ub, va, vb, err);
  if (err <= density * std::abs((ub - ua) * (vb - va)) || depth >= 12) return K;
  const double um = 0.5 * (ua + ub), vm = 0.5 * (va + vb);
  return integrate_2d(f, ua, um, va, vm, density, depth + 1) +
         integrate_2d(f, um, ub, va, vm, density, depth + 1) +
         integrate_2d(f, ua, um, vm, vb, density, depth + 1) +
         integrate_2d(f, um, ub, vm, vb, density, depth + 1);
}

// Break points of [a, b]: a, interior distinct knots, b.
std::vector<double> breaks(const KnotVector& U, double a, double b) {
  std::vector<double> x{a};
  for (double k : unique_knots(U).first)
    if (k > a + tol::kParametric && k < b - tol::kParametric) x.push_back(k);
  x.push_back(b);
  return x;
}

}  // namespace

// ------------------------------------------------------------------ curves

double curvature(const NurbsCurve& c, double u) {
  const auto D = c.derivatives(u, 2);
  const double n1 = norm(D[1]);
  if (n1 <= tol::kLength) throw GeometryError("curvature: C'(u) vanishes at u = " + std::to_string(u));
  return norm(cross(D[1], D[2])) / (n1 * n1 * n1);
}

double torsion(const NurbsCurve& c, double u) {
  const auto D = c.derivatives(u, 3);
  const Vec3 cr = cross(D[1], D[2]);
  const double n1 = norm(D[1]), ncr = norm(cr);
  if (n1 <= tol::kLength || ncr <= tol::kAngular * n1 * (n1 + norm(D[2])))
    throw GeometryError("torsion: undefined where the curvature vanishes (u = " + std::to_string(u) + ")");
  return dot(cr, D[3]) / (ncr * ncr);
}

double arc_length(const NurbsCurve& c, double u0, double u1, double rel_tol) {
  if (u0 == u1) return 0.0;
  if (u0 > u1) return -arc_length(c, u1, u0, rel_tol);
  const double a = c.domain_start(), b = c.domain_end();
  if (u0 < a - tol::kParametric || u1 > b + tol::kParametric)
    throw GeometryError("arc_length: parameters outside the curve domain");
  u0 = std::max(u0, a);
  u1 = std::min(u1, b);
  const auto speed = [&](double t) { return norm(c.derivatives(t, 1)[1]); };
  const auto x = breaks(c.knots(), u0, u1);
  // coarse estimate sets the absolute error budget
  double estimate = 0.0, err = 0.0;
  for (std::size_t i = 0; i + 1 < x.size(); ++i) estimate += gk15(speed, x[i], x[i + 1], err);
  const double density = rel_tol * std::max(estimate, tol::kLength) / (u1 - u0);
  double total = 0.0;
  for (std::size_t i = 0; i + 1 < x.size(); ++i) total += integrate_1d(speed, x[i], x[i + 1], density);
  return total;
}

double arc_length(const NurbsCurve& c, double rel_tol) {
  return arc_length(c, c.domain_start(), c.domain_end(), rel_tol);
}

double parameter_at_length(const NurbsCurve& c, double s, double rel_tol) {
  const double a = c.domain_start(), b = c.domain_end();
  const double L = arc_length(c, rel_tol);
  const double eps = rel_tol * std::max(L, tol::kLength);
  if (s < -eps || s > L + eps) throw GeometryError("parameter_at_length: s outside [0, length]");
  if (s <= 0.0) return a;
  if (s >= L) return b;
  // safeguarded Newton on F(u) = length(a, u) - s, F' = |C'(u)|; F tracked incrementally
  double lo = a, hi = b;
  double u = a + (b - a) * s / L;
  double F = arc_length(c, a, u, rel_tol) - s;
  for (int it = 0; it < 100 && std::abs(F) > eps; ++it) {
    if (F < 0.0) lo = u; else hi = u;
    const double speed = norm(c.derivatives(u, 1)[1]);
    double un = speed > 0.0 ? u - F / speed : 0.5 * (lo + hi);
    if (!(un > lo && un < hi)) un = 0.5 * (lo + hi);  // bisection fallback
    F += arc_length(c, u, un, rel_tol);
    u = un;
    if (hi - lo <= tol::kParametric) break;
  }
  return u;
}

// ------------------------------------------------------------------ surfaces

namespace {

// Scale-aware: a partial that is tiny relative to the other one (e.g. 1e-32 vs 5 at a pole,
// due to rounding) counts as vanishing, as do (nearly) parallel partials.
bool singular_point(const std::vector<std::vector<Vec3>>& d) {
  const double nu = norm(d[1][0]), nv = norm(d[0][1]);
  const double scale = std::max(nu, nv);
  return nu <= tol::kAngular * scale || nv <= tol::kAngular * scale ||
         norm(cross(d[1][0], d[0][1])) <= tol::kAngular * nu * nv;
}

// Curvature from the derivatives at a regular point.
SurfaceCurvature curvature_from(const std::vector<std::vector<Vec3>>& D) {
  const Vec3 &Su = D[1][0], &Sv = D[0][1];
  const Vec3 n = normalized(cross(Su, Sv));
  const double E = dot(Su, Su), F = dot(Su, Sv), G = dot(Sv, Sv);
  const double L = dot(D[2][0], n), M = dot(D[1][1], n), N = dot(D[0][2], n);
  const double det = E * G - F * F;

  SurfaceCurvature r;
  r.normal = n;
  r.gaussian = (L * N - M * M) / det;
  r.mean = (E * N - 2.0 * F * M + G * L) / (2.0 * det);
  const double disc = std::sqrt(std::max(r.mean * r.mean - r.gaussian, 0.0));
  r.k1 = r.mean + disc;
  r.k2 = r.mean - disc;

  // principal direction of k1: (L - kE) a + (M - kF) b = 0, (M - kF) a + (N - kG) b = 0
  const double p = L - r.k1 * E, q = M - r.k1 * F, t = N - r.k1 * G;
  // pick the better-conditioned equation; (a, b) = (-q, p) solves p a + q b = 0
  const Vec3 d = (std::abs(p) + std::abs(q) >= std::abs(q) + std::abs(t)) ? (-q) * Su + p * Sv
                                                                          : t * Su + (-q) * Sv;
  const double scale = std::abs(r.k1) + std::abs(r.k2) + 1e-300;
  // umbilic (every direction is principal) or numerically degenerate: use S_u
  const bool umbilic = disc <= 1e-8 * scale || norm(d) <= 1e-12 * (E + G) * scale;
  r.dir1 = umbilic ? normalized(Su) : normalized(d);
  r.dir2 = normalized(cross(n, r.dir1));
  return r;
}

}  // namespace

SurfaceCurvature curvature(const NurbsSurface& s, double u, double v) {
  const auto D = s.derivatives(u, v, 2);
  if (!singular_point(D)) return curvature_from(D);

  // Pole / collapsed edge: Richardson extrapolation X0 = 2 X(h) - X(2h) from two regular points
  // moved towards the interior (error O(h^2); rounding grows like eps / h^2, h = 1e-4 balances).
  const auto du = s.domain_u(), dv = s.domain_v();
  const auto at = [&](double f) {
    const double hu = f * (du[1] - du[0]), hv = f * (dv[1] - dv[0]);
    const double un = u < 0.5 * (du[0] + du[1]) ? u + hu : u - hu;
    const double vn = v < 0.5 * (dv[0] + dv[1]) ? v + hv : v - hv;
    const auto Dn = s.derivatives(un, vn, 2);
    if (singular_point(Dn))
      throw GeometryError("curvature: surface is singular near (" + std::to_string(u) + ", " +
                          std::to_string(v) + ")");
    return curvature_from(Dn);
  };
  const SurfaceCurvature c1 = at(1e-4), c2 = at(2e-4);
  SurfaceCurvature r;
  r.gaussian = 2.0 * c1.gaussian - c2.gaussian;
  r.mean = 2.0 * c1.mean - c2.mean;
  const double disc = std::sqrt(std::max(r.mean * r.mean - r.gaussian, 0.0));
  r.k1 = r.mean + disc;
  r.k2 = r.mean - disc;
  r.normal = s.normal(u, v);  // limit normal; same orientation as c1.normal
  const Vec3 d1 = c1.dir1 - dot(c1.dir1, r.normal) * r.normal;
  r.dir1 = normalized(d1);
  r.dir2 = normalized(cross(r.normal, r.dir1));
  return r;
}

double area(const NurbsSurface& s, double u0, double u1, double v0, double v1, double rel_tol) {
  const auto du = s.domain_u(), dv = s.domain_v();
  if (!(u0 < u1) || !(v0 < v1) || u0 < du[0] - tol::kParametric || u1 > du[1] + tol::kParametric ||
      v0 < dv[0] - tol::kParametric || v1 > dv[1] + tol::kParametric)
    throw GeometryError("area: parameter ranges must be increasing and inside the domain");
  u0 = std::max(u0, du[0]);
  u1 = std::min(u1, du[1]);
  v0 = std::max(v0, dv[0]);
  v1 = std::min(v1, dv[1]);
  const auto dA = [&](double a, double b) {
    const auto d = s.first_derivatives(a, b);
    return norm(cross(d[1], d[2]));
  };
  const auto xu = breaks(s.knots_u(), u0, u1), xv = breaks(s.knots_v(), v0, v1);
  double estimate = 0.0, err = 0.0;
  for (std::size_t i = 0; i + 1 < xu.size(); ++i)
    for (std::size_t j = 0; j + 1 < xv.size(); ++j)
      estimate += gk15_2d(dA, xu[i], xu[i + 1], xv[j], xv[j + 1], err);
  const double density = rel_tol * std::max(estimate, tol::kLength) / ((u1 - u0) * (v1 - v0));
  double total = 0.0;
  for (std::size_t i = 0; i + 1 < xu.size(); ++i)
    for (std::size_t j = 0; j + 1 < xv.size(); ++j)
      total += integrate_2d(dA, xu[i], xu[i + 1], xv[j], xv[j + 1], density);
  return total;
}

double area(const NurbsSurface& s, double rel_tol) {
  const auto du = s.domain_u(), dv = s.domain_v();
  return area(s, du[0], du[1], dv[0], dv[1], rel_tol);
}

}  // namespace nf
