#include "nurbforge/geometry/primitives.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>
#include <utility>
#include <vector>

#include "nurbforge/core/error.hpp"

namespace nf {

namespace detail {

std::size_t arc_segment_count(double sweep) {
  constexpr double half_pi = std::numbers::pi / 2.0;
  if (sweep <= half_pi + tol::kAngular) return 1;
  if (sweep <= std::numbers::pi + tol::kAngular) return 2;
  if (sweep <= 3.0 * half_pi + tol::kAngular) return 3;
  return 4;
}

KnotVector arc_knots(std::size_t narcs) {
  const std::size_t n = 2 * narcs;  // index of last control point
  KnotVector U(n + 4, 0.0);
  for (std::size_t i = n + 1; i < n + 4; ++i) U[i] = 1.0;
  for (std::size_t i = 1; i < narcs; ++i) {
    const double k = static_cast<double>(i) / static_cast<double>(narcs);
    U[1 + 2 * i] = k;
    U[2 + 2 * i] = k;
  }
  return U;
}

Vec3 intersect_lines(const Vec3& p0, const Vec3& t0, const Vec3& p2, const Vec3& t2) {
  // Solve p0 + a t0 = p2 + b t2 in the least-squares sense.
  const Vec3 w = p0 - p2;
  const double a = dot(t0, t0), b = dot(t0, t2), c = dot(t2, t2);
  const double d = dot(t0, w), e = dot(t2, w);
  const double den = a * c - b * b;
  if (std::abs(den) < tol::kLength) throw GeometryError("lines are parallel");
  const double s = (b * e - c * d) / den;
  return p0 + s * t0;
}

}  // namespace detail

Vec3 any_perpendicular(const Vec3& n) {
  const Vec3 u = normalized(n);
  const Vec3 helper = (std::abs(u.x) < 0.9) ? Vec3{1, 0, 0} : Vec3{0, 1, 0};
  return normalized(cross(u, helper));
}

NurbsCurve make_line(const Vec3& p0, const Vec3& p1) {
  return {1, {0.0, 0.0, 1.0, 1.0}, {p0, p1}};
}

NurbsCurve make_polyline(const std::vector<Vec3>& points) {
  if (points.size() < 2) throw GeometryError("polyline needs at least 2 points");
  // chord-length parametrisation keeps the speed roughly uniform
  KnotVector U(points.size() + 2);
  U[0] = 0.0;
  double total = 0.0;
  std::vector<double> acc(points.size(), 0.0);
  for (std::size_t i = 1; i < points.size(); ++i) {
    total += distance(points[i], points[i - 1]);
    acc[i] = total;
  }
  if (total < tol::kLength) throw GeometryError("polyline has zero length");
  for (std::size_t i = 0; i < points.size(); ++i) U[i + 1] = acc[i] / total;
  U.back() = 1.0;
  return {1, std::move(U), points};
}

NurbsCurve make_arc(const Vec3& O, const Vec3& x_axis, const Vec3& y_axis, double r,
                    double ths, double the) {
  if (r <= 0.0) throw GeometryError("radius must be positive");
  const Vec3 X = normalized(x_axis);
  const Vec3 Y = normalized(y_axis);
  if (std::abs(dot(X, Y)) > tol::kAngular) throw GeometryError("arc axes must be orthogonal");
  while (the < ths) the += 2.0 * std::numbers::pi;
  const double theta = the - ths;
  if (theta < tol::kAngular) throw GeometryError("arc sweep must be positive");

  const std::size_t narcs = detail::arc_segment_count(theta);
  const double dtheta = theta / static_cast<double>(narcs);
  const std::size_t n = 2 * narcs;
  const double w1 = std::cos(dtheta / 2.0);

  std::vector<Vec3> P(n + 1);
  std::vector<double> W(n + 1, 1.0);
  Vec3 P0 = O + r * std::cos(ths) * X + r * std::sin(ths) * Y;
  Vec3 T0 = -std::sin(ths) * X + std::cos(ths) * Y;
  P[0] = P0;
  std::size_t index = 0;
  double angle = ths;
  for (std::size_t i = 1; i <= narcs; ++i) {
    angle += dtheta;
    const Vec3 P2 = O + r * std::cos(angle) * X + r * std::sin(angle) * Y;
    const Vec3 T2 = -std::sin(angle) * X + std::cos(angle) * Y;
    P[index + 2] = P2;
    P[index + 1] = detail::intersect_lines(P0, T0, P2, T2);
    W[index + 1] = w1;
    index += 2;
    P0 = P2;
    T0 = T2;
  }
  return {2, detail::arc_knots(narcs), std::move(P), std::move(W)};
}

NurbsCurve make_circle(const Vec3& center, const Vec3& normal, double radius) {
  const Vec3 n = normalized(normal);
  const Vec3 X = any_perpendicular(n);
  const Vec3 Y = cross(n, X);
  return make_arc(center, X, Y, radius, 0.0, 2.0 * std::numbers::pi);
}

NurbsSurface make_bilinear(const Vec3& p00, const Vec3& p10, const Vec3& p01, const Vec3& p11) {
  return {1, 1, {0, 0, 1, 1}, {0, 0, 1, 1}, 2, 2, {p00, p01, p10, p11}};
}

NurbsSurface make_rectangle(const Vec3& o, const Vec3& eu, const Vec3& ev) {
  return make_bilinear(o, o + eu, o + ev, o + eu + ev);
}

// ------------------------------------------------------------------ constructors

namespace {

/// Centre and unit normal of the circle through three points (normal: (p1-p0) x (p2-p0)).
std::pair<Vec3, Vec3> circumcircle(const Vec3& p0, const Vec3& p1, const Vec3& p2) {
  const Vec3 a = p1 - p0, b = p2 - p0;
  const Vec3 n = cross(a, b);
  const double nn = dot(n, n);
  const double scale = std::max({dot(a, a), dot(b, b), 1e-300});
  if (nn <= 1e-24 * scale * scale)
    throw GeometryError("circle through 3 points: points are collinear");
  // circumcentre: p0 + (|a|^2 (b x n) + |b|^2 (n x a)) / (2 |n|^2)
  const Vec3 c = p0 + (1.0 / (2.0 * nn)) * (dot(a, a) * cross(b, n) + dot(b, b) * cross(n, a));
  return {c, (1.0 / std::sqrt(nn)) * n};
}

}  // namespace

NurbsCurve make_circle_3pt(const Vec3& p0, const Vec3& p1, const Vec3& p2) {
  const auto [c, n] = circumcircle(p0, p1, p2);
  const double r = distance(c, p0);
  const Vec3 x = (1.0 / r) * (p0 - c);
  return make_arc(c, x, cross(n, x), r, 0.0, 2.0 * std::numbers::pi);
}

NurbsCurve make_arc_3pt(const Vec3& p0, const Vec3& p1, const Vec3& p2) {
  const auto [c, n] = circumcircle(p0, p1, p2);
  const double r = distance(c, p0);
  const Vec3 x = (1.0 / r) * (p0 - c), y = cross(n, x);
  // p0, p1, p2 run counter-clockwise about n, so p2's angle is the sweep
  double end = std::atan2(dot(p2 - c, y), dot(p2 - c, x));
  if (end <= 0.0) end += 2.0 * std::numbers::pi;
  NurbsCurve arc = make_arc(c, x, y, r, 0.0, end);
  // snap the end points onto the given points (rounding of the trigonometry)
  auto P = arc.control_points();
  P.front() = p0;
  P.back() = p2;
  return {arc.degree(), arc.knots(), P, arc.weights()};
}

NurbsCurve make_polygon(const Vec3& center, const Vec3& normal, double radius, std::size_t n,
                        const Vec3& x_axis) {
  if (n < 3) throw GeometryError("polygon needs at least 3 vertices");
  if (!(radius > 0.0)) throw GeometryError("polygon radius must be > 0");
  if (!(norm(normal) > 0.0)) throw GeometryError("polygon normal must not be zero");
  const Vec3 z = normalized(normal);
  Vec3 x = x_axis - dot(x_axis, z) * z;
  x = norm(x) > 1e-12 * norm(x_axis) && norm(x) > 0.0 ? normalized(x) : any_perpendicular(z);
  const Vec3 y = cross(z, x);
  std::vector<Vec3> pts;
  for (std::size_t i = 0; i < n; ++i) {
    const double a = 2.0 * std::numbers::pi * static_cast<double>(i) / static_cast<double>(n);
    pts.push_back(center + radius * (std::cos(a) * x + std::sin(a) * y));
  }
  pts.push_back(pts.front());
  return make_polyline(pts);
}

NurbsCurve make_rectangle_curve(const Vec3& o, const Vec3& eu, const Vec3& ev) {
  if (!(norm(cross(eu, ev)) > 0.0)) throw GeometryError("rectangle edges must not be parallel");
  return make_polyline({o, o + eu, o + eu + ev, o + ev, o});
}

NurbsCurve make_helix(const Vec3& axis_point, const Vec3& axis, const Vec3& start, double pitch,
                      double turns, double end_radius, bool left_handed, double tolerance) {
  if (!(norm(axis) > 0.0)) throw GeometryError("helix: axis must not be zero");
  if (!(turns > 0.0) || !std::isfinite(turns)) throw GeometryError("helix: turns must be > 0");
  if (!(tolerance > 0.0)) throw GeometryError("helix: tolerance must be > 0");
  const Vec3 z = normalized(axis);
  const Vec3 rel = start - axis_point;
  const double h0 = dot(rel, z);
  const Vec3 radial = rel - h0 * z;
  const double r0 = norm(radial);
  const double r1 = end_radius < 0.0 ? r0 : end_radius;
  if (r0 == 0.0 && r1 == 0.0) throw GeometryError("helix: start and end radius are both 0");
  const Vec3 x = r0 > 1e-14 * (std::abs(h0) + 1.0) ? (1.0 / r0) * radial : any_perpendicular(z);
  const Vec3 y = left_handed ? cross(x, z) : cross(z, x);
  const double w = 2.0 * std::numbers::pi, dr = (r1 - r0) / turns;
  const Vec3 base = axis_point + h0 * z;
  // exact point and derivatives at t (turns)
  const auto eval = [&](double t, Vec3& p, Vec3& d1, Vec3& d2) {
    const double r = r0 + dr * t, c = std::cos(w * t), s = std::sin(w * t);
    const Vec3 e = c * x + s * y, f = -s * x + c * y;  // radial and tangential unit vectors
    p = base + (pitch * t) * z + r * e;
    d1 = pitch * z + dr * e + (r * w) * f;
    d2 = (2.0 * dr * w) * f - (r * w * w) * e;
  };
  const double scale = std::max({r0, r1, std::abs(pitch) * turns, 1e-300});
  if (tolerance < 1e-14 * scale) throw GeometryError("helix: tolerance below the rounding level");
  std::size_t n = std::max<std::size_t>(4, static_cast<std::size_t>(std::ceil(4.0 * turns)));
  for (int attempt = 0; attempt < 12 && n <= (std::size_t{1} << 22); ++attempt) {
    const double h = turns / static_cast<double>(n);
    KnotVector U(6, 0.0);
    std::vector<Vec3> P;
    P.reserve(5 * n + 1);
    Vec3 pa, da, ea;
    eval(0.0, pa, da, ea);
    P.push_back(pa);
    for (std::size_t i = 1; i <= n; ++i) {
      const double t1 = i == n ? turns : static_cast<double>(i) * h;
      Vec3 pb, db, eb;
      eval(t1, pb, db, eb);
      // quintic Hermite -> Bezier control points (h = span length in t)
      P.push_back(pa + (h / 5.0) * da);
      P.push_back(pa + (2.0 * h / 5.0) * da + (h * h / 20.0) * ea);
      P.push_back(pb - (2.0 * h / 5.0) * db + (h * h / 20.0) * eb);
      P.push_back(pb - (h / 5.0) * db);
      P.push_back(pb);
      for (int k = 0; k < (i == n ? 6 : 5); ++k) U.push_back(t1);
      pa = pb;
      da = db;
      ea = eb;
    }
    NurbsCurve c(5, U, P);
    // the Hermite data match up to the second derivative: the joints are C2, so two of the five
    // knot copies are removable (up to rounding)
    for (std::size_t i = 1; i < n; ++i)
      c = c.try_remove_knot(static_cast<double>(i) * h, 2, 1e-12 * scale).first;
    // verify against the exact curve at the same parameter
    double err = 0.0;
    for (std::size_t i = 0; i < n; ++i)
      for (int k = 1; k < 24; ++k) {
        const double t = (static_cast<double>(i) + k / 24.0) * h;
        Vec3 pe, d1, d2;
        eval(t, pe, d1, d2);
        err = std::max(err, distance(c.evaluate(t), pe));
      }
    if (err <= tolerance) return c;
    // the Hermite error falls like n^-6: aim at 0.75 tolerance (at least 1.25 n)
    const double grow = std::pow(err / (0.75 * tolerance), 1.0 / 6.0);
    n = std::max(n + n / 4 + 1, static_cast<std::size_t>(std::ceil(static_cast<double>(n) * grow)));
  }
  throw GeometryError("helix: tolerance not reached");
}

}  // namespace nf
