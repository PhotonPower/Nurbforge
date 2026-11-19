#include "nurbforge/modeling/blend.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

#include "nurbforge/core/error.hpp"
#include "nurbforge/geometry/primitives.hpp"
#include "nurbforge/modeling/construct.hpp"
#include "nurbforge/nurbs/knot_vector.hpp"

namespace nf {

namespace {

struct EndData {
  Vec3 p;  // end point
  Vec3 t;  // unit tangent, pointing along the blend direction
  Vec3 k;  // curvature vector (independent of the direction of travel)
};

EndData end_data(const NurbsCurve& c, CurveEnd e, bool leaving) {
  const double u = e == CurveEnd::End ? c.domain_end() : c.domain_start();
  const auto d = c.derivatives(u, 2);
  const double s = norm(d[1]);
  if (!(s > 0.0)) throw GeometryError("blend_curves: vanishing tangent at a curve end");
  EndData r;
  r.p = d[0];
  // a: leave the curve (forward at its end, backward at its start);
  // b: enter the curve towards its interior (forward at its start, backward at its end)
  const bool forward = leaving ? e == CurveEnd::End : e == CurveEnd::Start;
  r.t = (forward ? 1.0 : -1.0) / s * d[1];
  const Vec3 c2 = cross(d[1], d[2]);
  r.k = (1.0 / (s * s * s * s)) * cross(c2, d[1]);  // (C' x C'') x C' / |C'|^4
  return r;
}

}  // namespace

NurbsCurve blend_curves(const NurbsCurve& a, CurveEnd end0, const NurbsCurve& b, CurveEnd end1,
                        int continuity, double tension) {
  if (continuity != 1 && continuity != 2)
    throw GeometryError("blend_curves: continuity must be 1 (G1) or 2 (G2)");
  if (!(tension > 0.0)) throw GeometryError("blend_curves: tension must be positive");
  const EndData s = end_data(a, end0, true), e = end_data(b, end1, false);
  const double chord = distance(s.p, e.p);
  if (!(chord > 0.0)) throw GeometryError("blend_curves: the curve ends coincide");
  std::vector<Vec3> P;
  if (continuity == 1) {
    // cubic: C'(0) = 3 (P1 - P0) = m t0, m = chord
    const double m = tension * chord;
    P = {s.p, s.p + (m / 3.0) * s.t, e.p - (m / 3.0) * e.t, e.p};
  } else {
    // quintic: C'(0) = 5 (P1 - P0) = m t0, C''(0) = 20 (P2 - 2 P1 + P0) = m^2 k0 (G2 with
    // zero tangential acceleration); symmetric at the end
    const double m = tension * chord;
    const Vec3 p1 = s.p + (m / 5.0) * s.t, q1 = e.p - (m / 5.0) * e.t;
    const Vec3 p2 = 2.0 * p1 - s.p + (m * m / 20.0) * s.k;
    const Vec3 q2 = 2.0 * q1 - e.p + (m * m / 20.0) * e.k;
    P = {s.p, p1, p2, q2, q1, e.p};
  }
  const std::size_t p = P.size() - 1;
  return {p, clamped_uniform_knots(p, P.size()), P, {}};
}

std::vector<NurbsCurve> tween_curves(const NurbsCurve& a, const NurbsCurve& b, std::size_t count) {
  if (count == 0) throw GeometryError("tween_curves: count must be at least 1");
  std::vector<NurbsCurve> c{a, b};
  make_compatible(c);
  const auto& Pa = c[0].control_points();
  const auto& Pb = c[1].control_points();
  const bool rational = c[0].is_rational() || c[1].is_rational();
  const auto weight = [](const NurbsCurve& x, std::size_t i) {
    return x.weights().empty() ? 1.0 : x.weights()[i];
  };
  std::vector<NurbsCurve> out;
  for (std::size_t k = 1; k <= count; ++k) {
    const double s = static_cast<double>(k) / static_cast<double>(count + 1);
    std::vector<Vec3> P(Pa.size());
    std::vector<double> W;
    for (std::size_t i = 0; i < Pa.size(); ++i) {
      const double wa = weight(c[0], i), wb = weight(c[1], i);
      const double w = (1.0 - s) * wa + s * wb;
      P[i] = ((1.0 - s) * wa * Pa[i] + s * wb * Pb[i]) / w;
      if (rational) W.push_back(w);
    }
    out.emplace_back(c[0].degree(), c[0].knots(), std::move(P), std::move(W));
  }
  return out;
}

namespace {

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

double extent(const std::vector<const NurbsCurve*>& curves, const std::vector<Vec3>& extra) {
  Vec3 ref = extra.empty() ? curves.front()->control_points().front() : extra.front();
  double size = 0.0;
  for (const auto* c : curves)
    for (const auto& p : c->control_points()) size = std::max(size, distance(p, ref));
  for (const auto& p : extra) size = std::max(size, distance(p, ref));
  return std::max(size, 1e-300);
}

void check_planar(const std::vector<const NurbsCurve*>& curves, const std::vector<Vec3>& extra,
                  const Vec3& n, double size) {
  const Vec3 ref = extra.empty() ? curves.front()->control_points().front() : extra.front();
  const auto off = [&](const Vec3& p) { return std::abs(dot(p - ref, n)) > 1e-9 * size; };
  bool bad = std::any_of(extra.begin(), extra.end(), off);
  for (const auto* c : curves)
    bad = bad || std::any_of(c->control_points().begin(), c->control_points().end(), off);
  if (bad) throw GeometryError("tangent_lines: the input is not in a plane normal to plane_normal");
}

std::size_t samples(const NurbsCurve& c, std::size_t per_span, std::size_t lo, std::size_t hi) {
  const std::size_t spans = std::max<std::size_t>(1, unique_knots(c.knots()).first.size() - 1);
  return std::clamp(per_span * spans, lo, hi);
}

}  // namespace

std::vector<TangentLine> tangent_lines(const Vec3& point, const NurbsCurve& curve,
                                       const Vec3& plane_normal) {
  if (!(norm(plane_normal) > 0.0)) throw GeometryError("tangent_lines: plane normal is zero");
  const Vec3 n = normalized(plane_normal);
  const double size = extent({&curve}, {point});
  check_planar({&curve}, {point}, n, size);
  const auto f = [&](double t) {
    const auto D = curve.derivatives(t, 1);
    return dot(n, cross(D[0] - point, D[1]));
  };
  const double t0 = curve.domain_start(), t1 = curve.domain_end();
  const std::size_t N = samples(curve, 64, 64, 4096);
  std::vector<double> roots;
  double ta = t0, fa = f(t0);
  if (fa == 0.0) roots.push_back(t0);
  for (std::size_t i = 1; i <= N; ++i) {
    const double tb =
        i == N ? t1 : t0 + (t1 - t0) * static_cast<double>(i) / static_cast<double>(N);
    const double fb = f(tb);
    if (fb == 0.0) {
      roots.push_back(tb);
    } else if (fa != 0.0 && (fa < 0.0) != (fb < 0.0)) {
      double lo = ta, hi = tb, flo = fa;
      for (int it = 0; it < 200 && hi - lo > 1e-16 * (1.0 + std::abs(lo)); ++it) {
        const double m = 0.5 * (lo + hi), fm = f(m);
        if ((fm < 0.0) == (flo < 0.0)) {
          lo = m;
          flo = fm;
        } else {
          hi = m;
        }
      }
      roots.push_back(0.5 * (lo + hi));
    }
    ta = tb;
    fa = fb;
  }
  std::vector<TangentLine> out;
  for (const double t : roots) {
    const Vec3 q = curve.evaluate(t);
    if (distance(q, point) <= 1e-9 * size) continue;  // the point itself lies on the curve
    if (std::any_of(out.begin(), out.end(),
                    [&](const TangentLine& l) { return distance(l.point_b, q) <= 1e-9 * size; }))
      continue;  // duplicate (closed curve start/end, or a root on a sample point)
    out.push_back({make_line(point, q), point, q, kNaN, t});
  }
  return out;
}

std::vector<TangentLine> tangent_lines(const NurbsCurve& a, const NurbsCurve& b,
                                       const Vec3& plane_normal) {
  if (!(norm(plane_normal) > 0.0)) throw GeometryError("tangent_lines: plane normal is zero");
  const Vec3 n = normalized(plane_normal);
  const double size = extent({&a, &b}, {});
  check_planar({&a, &b}, {}, n, size);
  const double s0 = a.domain_start(), s1 = a.domain_end(), t0 = b.domain_start(),
               t1 = b.domain_end();
  const std::size_t Ns = samples(a, 32, 64, 400), Nt = samples(b, 32, 64, 400);
  std::vector<double> S(Ns + 1), T(Nt + 1);
  std::vector<Vec3> A(Ns + 1), dA(Ns + 1), B(Nt + 1), dB(Nt + 1);
  for (std::size_t i = 0; i <= Ns; ++i) {
    S[i] = i == Ns ? s1 : s0 + (s1 - s0) * static_cast<double>(i) / static_cast<double>(Ns);
    const auto D = a.derivatives(S[i], 1);
    A[i] = D[0];
    dA[i] = D[1];
  }
  for (std::size_t j = 0; j <= Nt; ++j) {
    T[j] = j == Nt ? t1 : t0 + (t1 - t0) * static_cast<double>(j) / static_cast<double>(Nt);
    const auto D = b.derivatives(T[j], 1);
    B[j] = D[0];
    dB[j] = D[1];
  }
  const auto F = [&](std::size_t i, std::size_t j) {
    const Vec3 d = B[j] - A[i];
    return std::array<double, 2>{dot(n, cross(d, dA[i])), dot(n, cross(d, dB[j]))};
  };
  std::vector<TangentLine> out;
  for (std::size_t i = 0; i < Ns; ++i)
    for (std::size_t j = 0; j < Nt; ++j) {
      const std::array<double, 2> c[4] = {F(i, j), F(i + 1, j), F(i, j + 1), F(i + 1, j + 1)};
      bool ok = true;
      for (std::size_t k = 0; k < 2 && ok; ++k) {
        const double mn = std::min({c[0][k], c[1][k], c[2][k], c[3][k]}),
                     mx = std::max({c[0][k], c[1][k], c[2][k], c[3][k]});
        ok = mn <= 0.0 && mx >= 0.0;
      }
      if (!ok) continue;
      // Newton from the cell centre
      double s = 0.5 * (S[i] + S[i + 1]), t = 0.5 * (T[j] + T[j + 1]);
      bool converged = false;
      for (int it = 0; it < 50; ++it) {
        const auto Da = a.derivatives(s, 2), Db = b.derivatives(t, 2);
        const Vec3 d = Db[0] - Da[0];
        const double f1 = dot(n, cross(d, Da[1])), f2 = dot(n, cross(d, Db[1]));
        const double scale = norm(d) + 1e-300;
        if (std::abs(f1) <= 1e-13 * scale * norm(Da[1]) &&
            std::abs(f2) <= 1e-13 * scale * norm(Db[1])) {
          converged = true;
          break;
        }
        const double j11 = dot(n, cross(d, Da[2])), j12 = dot(n, cross(Db[1], Da[1]));
        const double j21 = dot(n, cross(-1.0 * Da[1], Db[1])), j22 = dot(n, cross(d, Db[2]));
        const double det = j11 * j22 - j12 * j21;
        if (!(std::abs(det) > 1e-300)) break;
        const double ds = (-f1 * j22 + f2 * j12) / det, dt = (-f2 * j11 + f1 * j21) / det;
        s = std::clamp(s + ds, s0, s1);
        t = std::clamp(t + dt, t0, t1);
        if (std::abs(ds) <= 1e-15 * (s1 - s0) && std::abs(dt) <= 1e-15 * (t1 - t0)) {
          const auto Ea = a.derivatives(s, 1), Eb = b.derivatives(t, 1);
          const Vec3 e = Eb[0] - Ea[0];
          converged = std::abs(dot(n, cross(e, Ea[1]))) <= 1e-9 * norm(e) * norm(Ea[1]) &&
                      std::abs(dot(n, cross(e, Eb[1]))) <= 1e-9 * norm(e) * norm(Eb[1]);
          break;
        }
      }
      if (!converged || s < S[i] - (S[i + 1] - S[i]) || s > S[i + 1] + (S[i + 1] - S[i]))
        continue;  // converged far away: that root belongs to (and is found from) another cell
      const Vec3 p = a.evaluate(s), q = b.evaluate(t);
      if (distance(p, q) <= 1e-9 * size) continue;  // intersection point, not a tangent line
      const bool dup = std::any_of(out.begin(), out.end(), [&](const TangentLine& l) {
        return distance(l.point_a, p) <= 1e-8 * size && distance(l.point_b, q) <= 1e-8 * size;
      });
      if (!dup) out.push_back({make_line(p, q), p, q, s, t});
    }
  std::sort(out.begin(), out.end(), [](const TangentLine& x, const TangentLine& y) {
    return x.param_a != y.param_a ? x.param_a < y.param_a : x.param_b < y.param_b;
  });
  return out;
}

}  // namespace nf
