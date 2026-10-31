#include "nurbforge/geometry/primitives.hpp"

#include <cmath>
#include <numbers>

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

}  // namespace nf
