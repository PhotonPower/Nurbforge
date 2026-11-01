#include "nurbforge/modeling/offset.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <optional>
#include <string>
#include <vector>

#include "nurbforge/analysis/properties.hpp"
#include "nurbforge/core/error.hpp"
#include "nurbforge/fitting/fit.hpp"
#include "nurbforge/nurbs/join.hpp"

namespace nf {

namespace {

constexpr std::size_t kCurveDegree = 3;
constexpr std::size_t kSurfaceDegree = 5;
constexpr std::size_t kMaxBreaks = 2049;  // per curve segment / surface direction
constexpr double kCuspMargin = 1e-9;      // reject d * k >= 1 - margin

using Func1 = std::function<Vec3(double)>;
using Func2 = std::function<Vec3(double, double)>;

// Break points of a fitted spline and the knot multiplicity at each interior break.
struct Breaks {
  std::vector<double> t;
  std::vector<std::size_t> mult;
};

KnotVector knots_of(const Breaks& b, std::size_t p) {
  KnotVector U(p + 1, b.t.front());
  for (std::size_t i = 1; i + 1 < b.t.size(); ++i) U.insert(U.end(), b.mult[i], b.t[i]);
  U.insert(U.end(), p + 1, b.t.back());
  return U;
}

// Data parameters: p + 2 per span (shared span ends), enough for every multiplicity <= p.
std::vector<double> params_of(const Breaks& b, std::size_t p) {
  std::vector<double> t{b.t.front()};
  const std::size_t per = p + 2;
  for (std::size_t i = 0; i + 1 < b.t.size(); ++i)
    for (std::size_t q = 1; q <= per; ++q)
      t.push_back(q == per ? b.t[i + 1]
                           : b.t[i] + (b.t[i + 1] - b.t[i]) * static_cast<double>(q) /
                                          static_cast<double>(per));
  return t;
}

// Halves every flagged span (new breaks are simple knots).
void refine(Breaks& b, const std::vector<bool>& bad) {
  Breaks out{{b.t.front()}, {b.mult.front()}};
  for (std::size_t i = 0; i + 1 < b.t.size(); ++i) {
    if (bad[i]) {
      out.t.push_back(0.5 * (b.t[i] + b.t[i + 1]));
      out.mult.push_back(1);
    }
    out.t.push_back(b.t[i + 1]);
    out.mult.push_back(b.mult[i + 1]);
  }
  b = std::move(out);
}

double lerp(double a, double b, std::size_t q, std::size_t n) {
  return q == n ? b : a + (b - a) * static_cast<double>(q) / static_cast<double>(n);
}

std::string num(double x) { return std::to_string(x); }

double angle_between(const Vec3& a, const Vec3& b) {
  return std::atan2(norm(cross(a, b)), dot(a, b));
}

// Least-squares fit of f on the breaks; flags the spans where |fit - f| > tol at checks + 1
// points. Returns true if any span was flagged.
bool fit_and_flag(const Func1& f, const Breaks& br, std::size_t p, std::size_t checks, double tol,
                  std::vector<bool>& bad, std::optional<NurbsCurve>* out = nullptr) {
  const auto t = params_of(br, p);
  std::vector<Vec3> pts;
  pts.reserve(t.size());
  for (double x : t) pts.push_back(f(x));
  NurbsCurve c = approximate_curve(pts, p, knots_of(br, p), t);
  bool any = false;
  for (std::size_t i = 0; i + 1 < br.t.size(); ++i) {
    if (bad[i]) continue;
    for (std::size_t q = 0; q <= checks; ++q) {
      const double x = lerp(br.t[i], br.t[i + 1], q, checks);
      if (distance(c.evaluate(x), f(x)) > tol) {
        bad[i] = any = true;
        break;
      }
    }
  }
  if (out != nullptr) *out = std::move(c);
  return any;
}

// ------------------------------------------------------------------ curves

// Adaptive cubic fit of f on [a, b]: spans are halved where |fit - f| > tol at 17 check points.
NurbsCurve fit_piece(const Func1& f, double a, double b, double tol) {
  Breaks br{{a, b}, {1, 1}};
  while (true) {
    std::vector<bool> bad(br.t.size() - 1, false);
    std::optional<NurbsCurve> c;
    if (!fit_and_flag(f, br, kCurveDegree, 16, tol, bad, &c)) return std::move(*c);
    if (br.t.size() >= kMaxBreaks)
      throw GeometryError("offset_curve: tolerance " + num(tol) + " not reached on [" + num(a) +
                          ", " + num(b) + "] (offset close to a cusp?)");
    refine(br, bad);
  }
}

// Exact offset point of a Bézier segment, with the regularity checks.
Vec3 offset_point(const NurbsCurve& seg, double t, double d, const Vec3& n, double min_speed) {
  const auto D = seg.derivatives(t, 2);
  const double speed = norm(D[1]);
  if (!(speed > min_speed))
    throw GeometryError("offset_curve: the curve derivative vanishes at u = " + num(t) +
                        " (offset direction undefined)");
  if (d != 0.0) {
    const double kappa = curvature(seg, t);
    const double signed_kappa = dot(n, cross(D[1], D[2])) >= 0.0 ? kappa : -kappa;
    if (d * signed_kappa >= 1.0 - kCuspMargin)
      throw GeometryError("offset_curve: |distance| = " + num(std::abs(d)) +
                          " reaches the radius of curvature " + num(1.0 / kappa) + " at u = " +
                          num(t) + " on the side of the centre of curvature (the offset would "
                          "have a cusp)");
  }
  return D[0] + d * normalized(cross(n, D[1]));
}

// ------------------------------------------------------------------ surfaces

// Unique knot values of the domain.
std::vector<double> spans_of(const KnotVector& U, std::size_t p) {
  std::vector<double> out;
  const double a = U[p], b = U[U.size() - p - 1];
  for (double k : U)
    if (k >= a && k <= b && (out.empty() || k > out.back())) out.push_back(k);
  return out;
}

// Knot lines of s as initial break lines of the fit.
Breaks surface_breaks(const KnotVector& U, std::size_t p) {
  const auto [values, mults] = unique_knots(U);
  const double a = U[p], b = U[U.size() - p - 1];
  Breaks br;
  for (std::size_t i = 0; i < values.size(); ++i) {
    if (values[i] < a || values[i] > b) continue;
    // s is C^k across the knot (k = p - m, at least); its offset is C^(k-1), which the fit of
    // degree q reproduces with multiplicity q - (k - 1), clamped to [1, q].
    const int k = static_cast<int>(p) - mults[i];
    const int q = static_cast<int>(kSurfaceDegree);
    br.t.push_back(values[i]);
    br.mult.push_back(static_cast<std::size_t>(std::clamp(q + 1 - k, 1, q)));
  }
  return br;
}

// Normal jumps across interior knot lines (one-sided normals from the Bézier patches).
void check_creases(const NurbsSurface& s, double d, double tol) {
  if (d == 0.0) return;
  const auto xu = spans_of(s.knots_u(), s.degree_u()), xv = spans_of(s.knots_v(), s.degree_v());
  if (xu.size() <= 2 && xv.size() <= 2) return;
  const auto patches = s.decompose_bezier();
  const std::size_t nsu = xu.size() - 1, nsv = xv.size() - 1;
  constexpr std::size_t kN = 4;
  const auto check = [&](const NurbsSurface& p0, const NurbsSurface& p1, double u, double v) {
    if (std::abs(d) * angle_between(p0.normal(u, v), p1.normal(u, v)) > tol)
      throw GeometryError("offset_surface: the normal is discontinuous at (u, v) = (" + num(u) +
                          ", " + num(v) + "); the offset would have a gap");
  };
  for (std::size_t a = 0; a + 1 < nsu; ++a)
    for (std::size_t b = 0; b < nsv; ++b)
      for (std::size_t r = 0; r <= kN; ++r)
        check(patches[a * nsv + b], patches[(a + 1) * nsv + b], xu[a + 1],
              lerp(xv[b], xv[b + 1], r, kN));
  for (std::size_t a = 0; a < nsu; ++a)
    for (std::size_t b = 0; b + 1 < nsv; ++b)
      for (std::size_t q = 0; q <= kN; ++q)
        check(patches[a * nsv + b], patches[a * nsv + b + 1], lerp(xu[a], xu[a + 1], q, kN),
              xv[b + 1]);
}

// d * k >= 1 for a principal curvature k: the offset degenerates (focal surface reached).
void check_curvature(const NurbsSurface& s, double d) {
  if (d == 0.0) return;
  const auto xu = spans_of(s.knots_u(), s.degree_u()), xv = spans_of(s.knots_v(), s.degree_v());
  constexpr std::size_t kN = 8;
  for (std::size_t a = 0; a + 1 < xu.size(); ++a)
    for (std::size_t b = 0; b + 1 < xv.size(); ++b)
      for (std::size_t q = 0; q <= kN; ++q)
        for (std::size_t r = 0; r <= kN; ++r) {
          const double u = lerp(xu[a], xu[a + 1], q, kN), v = lerp(xv[b], xv[b + 1], r, kN);
          const SurfaceCurvature c = curvature(s, u, v);
          const double k = d > 0.0 ? c.k1 : c.k2;  // the one maximising d * k
          if (d * k >= 1.0 - kCuspMargin)
            throw GeometryError("offset_surface: |distance| = " + num(std::abs(d)) +
                                " reaches the radius of curvature " + num(1.0 / std::abs(k)) +
                                " at (u, v) = (" + num(u) + ", " + num(v) +
                                ") on the side of the centre of curvature (the offset would "
                                "degenerate; local self-intersections are not removed)");
        }
}

// 1D presizing of the breaks of one direction: fits along iso-lines at the given transversal
// parameters until every iso-line fit is within tol.
void presize(Breaks& br, const Func2& g, const std::vector<double>& across, double tol) {
  while (true) {
    std::vector<bool> bad(br.t.size() - 1, false);
    bool any = false;
    for (double w : across)
      any = fit_and_flag([&](double t) { return g(t, w); }, br, kSurfaceDegree, 10, tol, bad) || any;
    if (!any) return;
    if (br.t.size() >= kMaxBreaks)
      throw GeometryError("offset_surface: tolerance " + num(tol) +
                          " not reached (offset close to a focal point?)");
    refine(br, bad);
  }
}

std::vector<double> samples_of(const Breaks& br, std::size_t n) {
  std::vector<double> out;
  for (std::size_t i = 0; i + 1 < br.t.size(); ++i)
    for (std::size_t q = (i == 0 ? 0 : 1); q <= n; ++q) out.push_back(lerp(br.t[i], br.t[i + 1], q, n));
  return out;
}

// Adaptive tensor-product fit of f: break lines presized by iso-line fits, then a full fit
// (rows along u, then columns along v) verified at 9 x 9 points per knot rectangle;
// rectangles failing the check are refined in both directions.
NurbsSurface fit_tensor(const Func2& f, Breaks bu, Breaks bv, double tol) {
  constexpr std::size_t p = kSurfaceDegree, kChecks = 8;
  const auto su = samples_of(bu, 4), sv = samples_of(bv, 4);
  presize(bu, f, sv, 0.5 * tol);
  presize(bv, [&](double v, double u) { return f(u, v); }, su, 0.5 * tol);
  while (true) {
    const KnotVector U = knots_of(bu, p), V = knots_of(bv, p);
    const auto us = params_of(bu, p), vs = params_of(bv, p);
    std::vector<std::vector<Vec3>> R;  // R[j] = control points of the row v = vs[j]
    R.reserve(vs.size());
    for (double v : vs) {
      std::vector<Vec3> row;
      row.reserve(us.size());
      for (double u : us) row.push_back(f(u, v));
      R.push_back(approximate_curve(row, p, U, us).control_points());
    }
    const std::size_t nu = R.front().size();
    std::vector<Vec3> P;
    std::size_t nv = 0;
    for (std::size_t i = 0; i < nu; ++i) {
      std::vector<Vec3> col;
      col.reserve(R.size());
      for (const auto& r : R) col.push_back(r[i]);
      const auto c = approximate_curve(col, p, V, vs).control_points();
      nv = c.size();
      P.insert(P.end(), c.begin(), c.end());
    }
    NurbsSurface fit{p, p, U, V, nu, nv, std::move(P)};

    std::vector<bool> badu(bu.t.size() - 1, false), badv(bv.t.size() - 1, false);
    bool any = false;
    for (std::size_t a = 0; a + 1 < bu.t.size(); ++a)
      for (std::size_t b = 0; b + 1 < bv.t.size(); ++b) {
        bool bad = false;
        for (std::size_t q = 0; q <= kChecks && !bad; ++q)
          for (std::size_t r = 0; r <= kChecks && !bad; ++r) {
            const double u = lerp(bu.t[a], bu.t[a + 1], q, kChecks);
            const double v = lerp(bv.t[b], bv.t[b + 1], r, kChecks);
            bad = distance(fit.evaluate(u, v), f(u, v)) > tol;
          }
        if (bad) badu[a] = badv[b] = any = true;
      }
    if (!any) return fit;
    if (bu.t.size() >= kMaxBreaks || bv.t.size() >= kMaxBreaks)
      throw GeometryError("offset_surface: tolerance " + num(tol) +
                          " not reached (offset close to a focal point?)");
    refine(bu, badu);
    refine(bv, badv);
  }
}

void check_args(const char* who, double d, double tol) {
  if (!std::isfinite(d)) throw GeometryError(std::string(who) + ": distance must be finite");
  if (!(tol > 0.0)) throw GeometryError(std::string(who) + ": tolerance must be positive");
}

}  // namespace

NurbsCurve offset_curve(const NurbsCurve& c, double distance, const Vec3& plane_normal,
                        double tolerance) {
  check_args("offset_curve", distance, tolerance);
  const double nn = norm(plane_normal);
  if (!(nn > 0.0) || !std::isfinite(nn))
    throw GeometryError("offset_curve: plane normal must be a non-zero vector");
  const Vec3 n = plane_normal / nn;

  // A NURBS curve lies in a plane iff its control points do (the basis is linearly independent).
  const auto& P = c.control_points();
  double extent = 0.0;
  for (const Vec3& q : P) {
    const double h = std::abs(dot(q - P.front(), n));
    if (h > tolerance)
      throw GeometryError("offset_curve: the curve does not lie in a plane normal to the given "
                          "plane normal (control point off the plane by " + num(h) + ")");
    extent = std::max(extent, nf::distance(q, P.front()));
  }
  const double min_speed = 1e-12 * std::max(extent, 1e-300) / (c.domain_end() - c.domain_start());

  // One fit per Bézier segment (tolerance/2), gaps at the joints <= tolerance/2 (closed at
  // their midpoint by join), so the total deviation stays below the tolerance.
  const auto segs = c.decompose_bezier();
  std::optional<NurbsCurve> result;
  for (std::size_t k = 0; k < segs.size(); ++k) {
    const NurbsCurve& seg = segs[k];
    const double a = seg.domain_start(), b = seg.domain_end();
    if (k > 0) {
      const Vec3 t0 = segs[k - 1].derivatives(a, 1)[1], t1 = seg.derivatives(a, 1)[1];
      if (!(norm(t0) > min_speed) || !(norm(t1) > min_speed))
        throw GeometryError("offset_curve: the curve derivative vanishes at u = " + num(a) +
                            " (offset direction undefined)");
      if (std::abs(distance) * angle_between(t0, t1) > 0.5 * tolerance)
        throw GeometryError("offset_curve: tangent discontinuity (corner) at u = " + num(a) +
                            "; the offset would have a gap - split the curve there");
    }
    const auto f = [&](double t) { return offset_point(seg, t, distance, n, min_speed); };
    NurbsCurve piece = fit_piece(f, a, b, 0.5 * tolerance);
    result = result ? join(*result, piece, tolerance) : std::move(piece);
  }
  return *result;
}

NurbsSurface offset_surface(const NurbsSurface& s, double distance, double tolerance) {
  check_args("offset_surface", distance, tolerance);
  check_curvature(s, distance);
  check_creases(s, distance, tolerance);
  const auto f = [&](double u, double v) {
    const auto D = s.first_derivatives(u, v);
    const Vec3 n = cross(D[1], D[2]);
    const double scale = std::max(norm(D[1]), norm(D[2]));
    // regular point (same test as NurbsSurface::normal), otherwise the limit normal
    if (scale > 0.0 && norm(n) > tol::kAngular * scale * scale) return D[0] + (distance / norm(n)) * n;
    return D[0] + distance * s.normal(u, v);
  };
  return fit_tensor(f, surface_breaks(s.knots_u(), s.degree_u()),
                    surface_breaks(s.knots_v(), s.degree_v()), tolerance);
}

}  // namespace nf
