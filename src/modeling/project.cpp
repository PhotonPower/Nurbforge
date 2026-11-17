#include "nurbforge/modeling/project.hpp"

#include <algorithm>
#include <cmath>
#include <optional>
#include <vector>

#include "nurbforge/analysis/closest_point.hpp"
#include "nurbforge/analysis/intersect.hpp"
#include "nurbforge/core/error.hpp"
#include "nurbforge/fitting/fit.hpp"
#include "nurbforge/geometry/primitives.hpp"
#include "nurbforge/nurbs/knot_vector.hpp"

namespace nf {

namespace {

struct UV {
  double u, v;
};

// Point S(u, v) = P + s d by Newton from (u, v) (s eliminated by least squares), clamped to the
// domain. Returns nothing if it does not converge onto the line.
std::optional<UV> ray_newton(const NurbsSurface& S, const Vec3& P, const Vec3& d, UV w,
                             double tol) {
  const auto du = S.domain_u(), dv = S.domain_v();
  for (int it = 0; it < 50; ++it) {
    const auto D = S.first_derivatives(w.u, w.v);
    // residual perpendicular to the line: r = (X - P) - ((X - P) . d) d
    const Vec3 x = D[0] - P;
    const Vec3 r = x - dot(x, d) * d;
    if (norm(r) <= 0.1 * tol) return w;
    const Vec3 a = D[1] - dot(D[1], d) * d, b = D[2] - dot(D[2], d) * d;
    const double aa = dot(a, a), ab = dot(a, b), bb = dot(b, b), det = aa * bb - ab * ab;
    if (!(std::abs(det) > 1e-300)) return {};
    const double ra = dot(r, a), rb = dot(r, b);
    w.u = std::clamp(w.u - (bb * ra - ab * rb) / det, du[0], du[1]);
    w.v = std::clamp(w.v - (aa * rb - ab * ra) / det, dv[0], dv[1]);
  }
  const Vec3 x = S.evaluate(w.u, w.v) - P;
  if (norm(x - dot(x, d) * d) <= tol) return w;
  return {};
}

}  // namespace

ProjectedCurve project_curve(const NurbsCurve& curve, const NurbsSurface& S, const Vec3& dir,
                             double tolerance) {
  if (!(tolerance > 0.0)) throw GeometryError("project_curve: tolerance must be positive");
  const bool along = norm(dir) > 0.0;
  const Vec3 d = along ? normalized(dir) : Vec3{};
  const auto du = S.domain_u(), dv = S.domain_v();
  // model size for the rays
  double L = 0.0;
  for (const auto& p : S.control_points()) L = std::max(L, distance(p, S.control_points()[0]));
  for (const auto& p : curve.control_points()) L = std::max(L, distance(p, S.control_points()[0]));
  L = 2.0 * L + 1.0;

  std::optional<UV> prev;
  const auto project = [&](double t) -> UV {
    const Vec3 P = curve.evaluate(t);
    if (!along) {
      SurfacePoint sp = closest_point(S, P);
      if (prev) {  // prefer the branch next to the previous point
        const SurfacePoint loc = closest_point_local(S, P, prev->u, prev->v);
        if (loc.distance <= sp.distance + tolerance) sp = loc;
      }
      return {sp.u, sp.v};
    }
    if (prev)
      if (auto w = ray_newton(S, P, d, *prev, tolerance)) return *w;
    const auto line = make_line(P - L * d, P + L * d);
    const auto hits = intersect(line, S);
    if (hits.empty()) throw GeometryError("project_curve: the projection ray misses the surface");
    const auto best = std::min_element(hits.begin(), hits.end(), [&](const auto& x, const auto& y) {
      return std::abs(x.t - 0.5) < std::abs(y.t - 0.5);
    });
    return {best->u, best->v};
  };

  // initial parameters: 8 per knot span of the input curve
  std::vector<double> T;
  const auto k = unique_knots(curve.knots()).first;
  for (std::size_t i = 0; i + 1 < k.size(); ++i)
    for (int j = 0; j < 8; ++j) T.push_back(k[i] + (k[i + 1] - k[i]) * j / 8.0);
  T.push_back(curve.domain_end());

  for (int round = 0; round < 12; ++round) {
    std::vector<Vec3> W, X;
    prev.reset();
    for (const double t : T) {
      const UV w = project(t);
      if (prev && (std::abs(w.u - prev->u) > 0.5 * (du[1] - du[0]) ||
                   std::abs(w.v - prev->v) > 0.5 * (dv[1] - dv[0])))
        throw GeometryError(
            "project_curve: the projection crosses the seam of the surface; split "
            "the curve there");
      prev = w;
      W.push_back({w.u, w.v, 0.0});
      X.push_back(S.evaluate(w.u, w.v));
    }
    const std::size_t deg = std::min<std::size_t>(3, T.size() - 1);
    ProjectedCurve r{interpolate_curve(X, deg, T), interpolate_curve(W, deg, T)};
    // verify at the midpoints; refine where the error is too large
    std::vector<double> add;
    prev.reset();
    for (std::size_t i = 0; i + 1 < T.size(); ++i) {
      const double tm = 0.5 * (T[i] + T[i + 1]);
      prev = UV{W[i].x, W[i].y};
      const UV w = project(tm);
      const Vec3 exact = S.evaluate(w.u, w.v);
      const Vec3 pc = r.pcurve.evaluate(tm);
      const double e =
          std::max(distance(r.curve.evaluate(tm), exact), distance(S.evaluate(pc.x, pc.y), exact));
      if (e > tolerance) add.push_back(tm);
    }
    if (add.empty()) return r;
    if (T.size() + add.size() > 8193)
      throw GeometryError("project_curve: tolerance not reached (projection not smooth?)");
    T.insert(T.end(), add.begin(), add.end());
    std::sort(T.begin(), T.end());
  }
  throw GeometryError("project_curve: tolerance not reached");
}

}  // namespace nf
