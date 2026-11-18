#include "nurbforge/modeling/project.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <map>
#include <optional>
#include <string>
#include <vector>

#include "geodesic.hpp"
#include "nurbforge/analysis/closest_point.hpp"
#include "nurbforge/analysis/intersect.hpp"
#include "nurbforge/core/error.hpp"
#include "nurbforge/fitting/fit.hpp"
#include "nurbforge/geometry/primitives.hpp"
#include "nurbforge/geometry/transform.hpp"
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

using ExactUV = std::function<UV(double t, const std::optional<UV>& prev)>;

// Curve on a surface through exactly computed points: interpolates the parameter points
// (p-curve) and their images (3D curve) at the parameters T (cubic) and adds the midpoints
// where either deviates from the exact point by more than the tolerance.
ProjectedCurve fit_on_surface(const NurbsSurface& S, std::vector<double> T, const ExactUV& exact,
                              double tolerance, const std::string& who) {
  const auto du = S.domain_u(), dv = S.domain_v();
  for (int round = 0; round < 12; ++round) {
    std::vector<Vec3> W, X;
    std::optional<UV> prev;
    for (const double t : T) {
      const UV w = exact(t, prev);
      if (prev && (std::abs(w.u - prev->u) > 0.5 * (du[1] - du[0]) ||
                   std::abs(w.v - prev->v) > 0.5 * (dv[1] - dv[0])))
        throw GeometryError(who +
                            ": the curve crosses the seam of the surface; split the curve there");
      prev = w;
      W.push_back({w.u, w.v, 0.0});
      X.push_back(S.evaluate(w.u, w.v));
    }
    const std::size_t deg = std::min<std::size_t>(3, T.size() - 1);
    ProjectedCurve r{interpolate_curve(X, deg, T), interpolate_curve(W, deg, T)};
    std::vector<double> add;
    for (std::size_t i = 0; i + 1 < T.size(); ++i) {
      const double tm = 0.5 * (T[i] + T[i + 1]);
      const UV w = exact(tm, UV{W[i].x, W[i].y});
      const Vec3 ex = S.evaluate(w.u, w.v);
      const Vec3 pc = r.pcurve.evaluate(tm);
      const double e =
          std::max(distance(r.curve.evaluate(tm), ex), distance(S.evaluate(pc.x, pc.y), ex));
      if (e > tolerance) add.push_back(tm);
    }
    if (add.empty()) return r;
    if (T.size() + add.size() > 8193)
      throw GeometryError(who + ": tolerance not reached (curve not smooth?)");
    T.insert(T.end(), add.begin(), add.end());
    std::sort(T.begin(), T.end());
  }
  throw GeometryError(who + ": tolerance not reached");
}

// 8 parameters per knot span of the curve, plus the end.
std::vector<double> span_samples(const NurbsCurve& c) {
  std::vector<double> T;
  const auto k = unique_knots(c.knots()).first;
  for (std::size_t i = 0; i + 1 < k.size(); ++i)
    for (int j = 0; j < 8; ++j) T.push_back(k[i] + (k[i + 1] - k[i]) * j / 8.0);
  T.push_back(c.domain_end());
  return T;
}

// (a, b) with a S_u + b S_v = w (least squares) at (u, v).
UV lift(const std::array<Vec3, 3>& D, const Vec3& w) {
  const double aa = dot(D[1], D[1]), ab = dot(D[1], D[2]), bb = dot(D[2], D[2]);
  const double det = aa * bb - ab * ab;
  if (!(std::abs(det) > 1e-300 * std::max(1.0, aa * bb)))
    throw GeometryError("geodesic: degenerate surface (S_u x S_v = 0)");
  const double ra = dot(D[1], w), rb = dot(D[2], w);
  return {(bb * ra - ab * rb) / det, (aa * rb - ab * ra) / det};
}

// Geodesic ODE in the parameter domain, state (u, v, u', v') for arc length: the tangential
// part of S_uu u'^2 + 2 S_uv u'v' + S_vv v'^2 + S_u u'' + S_v v'' vanishes.
std::array<double, 4> geodesic_rhs(const NurbsSurface& S, const std::array<double, 4>& y) {
  const auto du = S.domain_u(), dv = S.domain_v();
  if (y[0] < du[0] - 1e-12 || y[0] > du[1] + 1e-12 || y[1] < dv[0] - 1e-12 || y[1] > dv[1] + 1e-12)
    throw GeometryError("the geodesic (offset) leaves the surface");
  const auto D = S.derivatives(std::clamp(y[0], du[0], du[1]), std::clamp(y[1], dv[0], dv[1]), 2);
  const Vec3 q = y[2] * y[2] * D[2][0] + 2.0 * y[2] * y[3] * D[1][1] + y[3] * y[3] * D[0][2];
  const UV acc = lift({D[0][0], D[1][0], D[0][1]}, -1.0 * q);
  return {y[2], y[3], acc.u, acc.v};
}

UV geodesic_end(const NurbsSurface& S, UV start, UV velocity, double length, int steps) {
  std::array<double, 4> y{start.u, start.v, velocity.u, velocity.v};
  const double h = length / steps;
  const auto axpy = [](const std::array<double, 4>& a, double f, const std::array<double, 4>& b) {
    return std::array<double, 4>{a[0] + f * b[0], a[1] + f * b[1], a[2] + f * b[2],
                                 a[3] + f * b[3]};
  };
  for (int i = 0; i < steps; ++i) {  // classical Runge-Kutta
    const auto k1 = geodesic_rhs(S, y);
    const auto k2 = geodesic_rhs(S, axpy(y, 0.5 * h, k1));
    const auto k3 = geodesic_rhs(S, axpy(y, 0.5 * h, k2));
    const auto k4 = geodesic_rhs(S, axpy(y, h, k3));
    for (std::size_t j = 0; j < 4; ++j) y[j] += h / 6.0 * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j]);
  }
  return {y[0], y[1]};
}

}  // namespace

ProjectedCurve project_curve(const NurbsCurve& curve, const NurbsSurface& S, const Vec3& dir,
                             double tolerance) {
  if (!(tolerance > 0.0)) throw GeometryError("project_curve: tolerance must be positive");
  const bool along = norm(dir) > 0.0;
  const Vec3 d = along ? normalized(dir) : Vec3{};
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

  return fit_on_surface(
      S, span_samples(curve),
      [&](double t, const std::optional<UV>& p) {
        prev = p;
        return project(t);
      },
      tolerance, "project_curve");
}

namespace {

Transform projection(const Plane& plane, const Vec3& direction) {
  const double nn = norm(plane.normal);
  if (!(nn > 0.0)) throw GeometryError("project_to_plane: the plane normal is zero");
  const Vec3 n = plane.normal / nn;
  const Vec3 d = norm(direction) > 0.0 ? normalized(direction) : n;
  const double dn = dot(d, n);
  if (std::abs(dn) <= tol::kAngular)
    throw GeometryError("project_to_plane: the direction is parallel to the plane");
  // x' = x - d (n.(x - p)) / (d.n) = (I - d n^T / dn) x + d (n.p) / dn
  const double c = dot(n, plane.point) / dn;
  return Transform({1 - d.x * n.x / dn, -d.x * n.y / dn, -d.x * n.z / dn, d.x * c, -d.y * n.x / dn,
                    1 - d.y * n.y / dn, -d.y * n.z / dn, d.y * c, -d.z * n.x / dn, -d.z * n.y / dn,
                    1 - d.z * n.z / dn, d.z * c});
}

}  // namespace

NurbsCurve project_to_plane(const NurbsCurve& curve, const Plane& plane, const Vec3& direction) {
  return curve.transformed(projection(plane, direction));
}

NurbsSurface project_to_plane(const NurbsSurface& surface, const Plane& plane,
                              const Vec3& direction) {
  return surface.transformed(projection(plane, direction));
}

ProjectedCurve interpolate_on_surface(const NurbsSurface& S, const std::vector<Vec3>& points,
                                      std::size_t degree, bool points_are_uv, double tolerance) {
  const std::string who = "interpolate_on_surface";
  if (points.size() < 2) throw GeometryError(who + ": needs at least two points");
  if (!(tolerance > 0.0)) throw GeometryError(who + ": tolerance must be positive");
  const auto du = S.domain_u(), dv = S.domain_v();
  std::vector<Vec3> uv, xyz;
  std::optional<UV> prev;
  for (const auto& p : points) {
    UV w{p.x, p.y};
    if (points_are_uv) {
      if (p.x < du[0] - 1e-12 || p.x > du[1] + 1e-12 || p.y < dv[0] - 1e-12 || p.y > dv[1] + 1e-12)
        throw GeometryError(who + ": a parameter point lies outside the surface domain");
    } else {
      SurfacePoint sp = closest_point(S, p);
      if (prev) {  // stay on the branch of the previous point
        const SurfacePoint loc = closest_point_local(S, p, prev->u, prev->v);
        if (loc.distance <= sp.distance + tolerance) sp = loc;
      }
      w = {sp.u, sp.v};
    }
    if (prev && (std::abs(w.u - prev->u) > 0.5 * (du[1] - du[0]) ||
                 std::abs(w.v - prev->v) > 0.5 * (dv[1] - dv[0])))
      throw GeometryError(who + ": consecutive points lie across the seam of the surface");
    prev = w;
    uv.push_back({w.u, w.v, 0.0});
    xyz.push_back(S.evaluate(w.u, w.v));
  }
  const std::size_t deg = std::min(degree, points.size() - 1);
  const NurbsCurve pc =
      interpolate_curve(uv, deg, data_parameters(xyz, Parametrization::ChordLength));
  ProjectedCurve r = fit_on_surface(
      S, span_samples(pc),
      [&](double t, const std::optional<UV>&) {
        const Vec3 q = pc.evaluate(t);
        return UV{std::clamp(q.x, du[0], du[1]), std::clamp(q.y, dv[0], dv[1])};
      },
      tolerance, who);
  r.pcurve = pc;  // the interpolating p-curve itself defines the curve
  return r;
}

ProjectedCurve offset_on_surface(const NurbsSurface& S, const NurbsCurve& pcurve, double dist,
                                 double tolerance) {
  const std::string who = "offset_on_surface";
  if (!(tolerance > 0.0)) throw GeometryError(who + ": tolerance must be positive");
  if (dist == 0.0) throw GeometryError(who + ": the distance must be non-zero");
  int last_steps = 8;
  const auto exact = [&](double t, const std::optional<UV>&) -> UV {
    const auto P = pcurve.derivatives(t, 1);
    const auto D = S.first_derivatives(P[0].x, P[0].y);
    const Vec3 T = P[1].x * D[1] + P[1].y * D[2];
    const Vec3 N = cross(D[1], D[2]);
    if (!(norm(T) > 0.0) || !(norm(N) > 0.0))
      throw GeometryError(who + ": degenerate curve or surface at t = " + std::to_string(t));
    const Vec3 B = normalized(cross(N, T)) * (dist > 0.0 ? 1.0 : -1.0);
    const auto w =
        detail::geodesic_point(S, P[0].x, P[0].y, B, std::abs(dist), tolerance, last_steps);
    return UV{w[0], w[1]};
  };
  std::map<double, UV> cache;  // fit_on_surface revisits the parameters of earlier rounds
  const auto cached = [&](double t, const std::optional<UV>& prev) {
    const auto it = cache.find(t);
    if (it != cache.end()) return it->second;
    return cache[t] = exact(t, prev);
  };
  return fit_on_surface(S, span_samples(pcurve), cached, tolerance, who);
}

namespace detail {

std::array<double, 2> geodesic_point(const NurbsSurface& S, double u, double v, const Vec3& dir,
                                     double length, double tolerance, int& steps) {
  const auto D = S.first_derivatives(u, v);
  const Vec3 N = cross(D[1], D[2]);
  if (!(norm(N) > 0.0)) throw GeometryError("geodesic: degenerate surface point");
  const Vec3 n = normalized(N);
  Vec3 d = dir - dot(dir, n) * n;
  if (!(norm(d) > 0.0)) throw GeometryError("geodesic: direction normal to the surface");
  d = normalized(d);
  const UV start{u, v}, vel = lift(D, d);
  int k = std::max(2, steps / 2);
  UV w = geodesic_end(S, start, vel, length, k);
  for (; k <= 4096; k *= 2) {
    const UV w2 = geodesic_end(S, start, vel, length, 2 * k);
    const double change = distance(S.evaluate(w.u, w.v), S.evaluate(w2.u, w2.v));
    w = w2;
    if (change <= 0.05 * tolerance) {
      steps = 2 * k;
      return {w.u, w.v};
    }
  }
  throw GeometryError("geodesic: no convergence");
}

}  // namespace detail

}  // namespace nf
