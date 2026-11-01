#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <algorithm>
#include <cmath>
#include <numbers>

#include "nurbforge/analysis/closest_point.hpp"
#include "nurbforge/core/error.hpp"
#include "nurbforge/fitting/fit.hpp"
#include "nurbforge/geometry/primitives.hpp"
#include "nurbforge/modeling/construct.hpp"
#include "nurbforge/modeling/offset.hpp"

using namespace nf;
constexpr double pi = std::numbers::pi;

namespace {

constexpr Vec3 ez{0, 0, 1};

// max |o(t) - c(t) - d * normalize(n x c'(t))| on a dense grid (parameter correspondence)
double curve_deviation(const NurbsCurve& c, const NurbsCurve& o, double d, const Vec3& n) {
  double e = 0.0;
  for (int i = 0; i <= 1000; ++i) {
    const double t = c.domain_start() + (c.domain_end() - c.domain_start()) * i / 1000.0;
    const Vec3 ref = c.evaluate(t) + d * normalized(cross(n, c.derivatives(t, 1)[1]));
    e = std::max(e, distance(o.evaluate(t), ref));
  }
  return e;
}

double surface_deviation(const NurbsSurface& s, const NurbsSurface& o, double d) {
  double e = 0.0;
  const auto du = s.domain_u(), dv = s.domain_v();
  for (int i = 0; i <= 60; ++i)
    for (int j = 0; j <= 60; ++j) {
      const double u = du[0] + (du[1] - du[0]) * i / 60.0, v = dv[0] + (dv[1] - dv[0]) * j / 60.0;
      e = std::max(e, distance(o.evaluate(u, v), s.evaluate(u, v) + d * s.normal(u, v)));
    }
  return e;
}

// max | dist(o(u,v), centre) - r | where centre(P) is computed from the point itself
template <class F>
double radial_error(const NurbsSurface& o, F&& signed_radius_error) {
  double e = 0.0;
  const auto du = o.domain_u(), dv = o.domain_v();
  for (int i = 0; i <= 60; ++i)
    for (int j = 0; j <= 60; ++j) {
      const double u = du[0] + (du[1] - du[0]) * i / 60.0, v = dv[0] + (dv[1] - dv[0]) * j / 60.0;
      e = std::max(e, std::abs(signed_radius_error(o.evaluate(u, v))));
    }
  return e;
}

NurbsSurface sphere(double R) {
  return revolve(make_arc({0, 0, 0}, {1, 0, 0}, {0, 0, 1}, R, -pi / 2, pi / 2), {0, 0, 0}, ez, 2 * pi);
}

// +1 if the surface normal at (u, v) points away from `centre`
double orientation(const NurbsSurface& s, double u, double v, const Vec3& centre) {
  return dot(s.normal(u, v), s.evaluate(u, v) - centre) > 0.0 ? 1.0 : -1.0;
}

}  // namespace

TEST_CASE("offset of a circle is a concentric circle", "[offset]") {
  const double R = 2.0;
  const Vec3 C{1, -1, 0.5};
  const auto c = make_circle(C, ez, R);
  const double d = GENERATE(0.5, -0.7, 1.5);
  const auto o = offset_curve(c, d, ez);
  CHECK(o.degree() == 3);
  CHECK(!o.is_rational());
  CHECK(o.domain_start() == c.domain_start());
  CHECK(o.domain_end() == c.domain_end());
  CHECK(curve_deviation(c, o, d, ez) <= 1e-7);
  // n x c' points inwards for a counter-clockwise circle about +z: radius R - d
  const double ccw = dot(cross(c.evaluate(0.0) - C, c.derivatives(0.0, 1)[1]), ez) > 0 ? 1.0 : -1.0;
  for (int i = 0; i <= 500; ++i) {
    const double t = c.domain_start() + (c.domain_end() - c.domain_start()) * i / 500.0;
    const Vec3 p = o.evaluate(t);
    CHECK(std::abs(distance(p, C) - (R - ccw * d)) <= 1e-7);
    CHECK(std::abs(p.z - C.z) <= 1e-12);
    // same polar angle as c(t)
    CHECK(norm(cross(normalized(p - C), normalized(c.evaluate(t) - C))) <= 1e-7);
  }
  CHECK(distance(o.evaluate(o.domain_start()), o.evaluate(o.domain_end())) <= 1e-12);  // closed
}

TEST_CASE("offset of a circular arc and tighter tolerances", "[offset]") {
  const double R = 3.0;
  const auto arc = make_arc({0, 0, 0}, {1, 0, 0}, {0, 1, 0}, R, 0.2, 2.5);
  for (double tol : {1e-5, 1e-9}) {
    const auto o = offset_curve(arc, -1.0, ez, tol);  // outwards
    CHECK(curve_deviation(arc, o, -1.0, ez) <= tol);
    for (int i = 0; i <= 300; ++i) {
      const double t = arc.domain_start() + (arc.domain_end() - arc.domain_start()) * i / 300.0;
      CHECK(std::abs(norm(o.evaluate(t)) - (R + 1.0)) <= tol);
    }
  }
  // the same arc seen from -z: offsets to the other side
  const auto o2 = offset_curve(arc, -1.0, {0, 0, -1});
  CHECK(std::abs(norm(o2.evaluate(1.0)) - (R - 1.0)) <= 1e-7);
}

TEST_CASE("offset of a line is a parallel line", "[offset]") {
  const auto l = make_line({0, 0, 0}, {4, 3, 0});
  const auto o = offset_curve(l, 2.0, ez);
  CHECK(distance(o.evaluate(o.domain_start()), {-1.2, 1.6, 0}) <= 1e-12);
  CHECK(distance(o.evaluate(o.domain_end()), {2.8, 4.6, 0}) <= 1e-12);
  CHECK(distance(o.evaluate(0.5), {0.8, 3.1, 0}) <= 1e-12);
  CHECK(o.num_control_points() == 4);  // a single cubic span reproduces it exactly
  // a polyline with a corner is rejected (the offset would have a gap)
  CHECK_THROWS_AS(offset_curve(make_polyline({{0, 0, 0}, {1, 0, 0}, {1, 1, 0}}), 0.1, ez),
                  GeometryError);
}

TEST_CASE("offset of a planar cubic matches the brute-force distance", "[offset]") {
  // tilted plane: the curve lies in the plane through the origin with normal n
  const Vec3 n = normalized({0.3, -0.2, 1.0});
  const Vec3 e1 = normalized(cross(n, {1, 0, 0})), e2 = cross(n, e1);
  std::vector<Vec3> pts;
  for (double x : {0.0, 1.0, 2.0, 3.0, 4.0, 5.0})
    pts.push_back(x * e1 + (0.6 * std::sin(1.1 * x)) * e2);
  const auto c = interpolate_curve(pts, 3);
  const double d = 0.25;
  const auto o = offset_curve(c, d, n);
  CHECK(curve_deviation(c, o, d, n) <= 1e-7);
  for (int i = 0; i <= 100; ++i) {
    const double t = c.domain_start() + (c.domain_end() - c.domain_start()) * i / 100.0;
    const Vec3 p = o.evaluate(t);
    CHECK(std::abs(dot(p, n)) <= 1e-9);                                    // stays in the plane
    CHECK(std::abs(closest_point(c, p).distance - d) <= 2e-7);             // at distance |d|
  }
}

TEST_CASE("offset_curve rejects singular and invalid input", "[offset]") {
  const auto c = make_circle({0, 0, 0}, ez, 1.0);
  CHECK_NOTHROW(offset_curve(c, 0.9, ez, 1e-6));
  CHECK_THROWS_AS(offset_curve(c, 1.0, ez), GeometryError);  // offset collapses to the centre
  CHECK_THROWS_AS(offset_curve(c, 1.5, ez), GeometryError);  // cusps / inverted
  CHECK_NOTHROW(offset_curve(c, -5.0, ez));                  // convex side: any distance
  CHECK_THROWS_AS(offset_curve(c, 0.1, {0, 0, 0}), GeometryError);
  CHECK_THROWS_AS(offset_curve(c, 0.1, {1, 0, 0}), GeometryError);  // not in that plane
  const auto spatial = NurbsCurve::from_points(3, {{0, 0, 0}, {1, 1, 0}, {2, 0, 1}, {3, 1, 0}});
  CHECK_THROWS_AS(offset_curve(spatial, 0.1, ez), GeometryError);
  // vanishing derivative (doubled control point at the start of a degree-1 curve)
  const NurbsCurve stall{1, {0, 0, 0.5, 1, 1}, {{0, 0, 0}, {0, 0, 0}, {1, 0, 0}}};
  CHECK_THROWS_AS(offset_curve(stall, 0.1, ez), GeometryError);
}

TEST_CASE("offset of a sphere is a concentric sphere", "[offset]") {
  const double R = 2.0, tol = 1e-6;
  const auto s = sphere(R);
  const double sgn = orientation(s, 0.3, 0.4, {0, 0, 0});
  for (double d : {0.5, -0.8}) {
    const auto o = offset_surface(s, d, tol);
    CHECK(o.degree_u() == 5);
    CHECK(o.degree_v() == 5);
    CHECK(!o.knots_u().empty());
    CHECK(o.domain_u() == s.domain_u());
    CHECK(o.domain_v() == s.domain_v());
    CHECK(surface_deviation(s, o, d) <= tol);
    CHECK(radial_error(o, [&](const Vec3& p) { return norm(p) - (R + sgn * d); }) <= tol);
    // the poles stay poles: the offset of the pole is the pole of the offset sphere
    CHECK(distance(o.evaluate(0.1, s.domain_v()[1]), {0, 0, R + sgn * d}) <= tol);
  }
  // offsetting into the centre (or beyond) is rejected
  CHECK_THROWS_AS(offset_surface(s, -sgn * R), GeometryError);
  CHECK_THROWS_AS(offset_surface(s, -sgn * 3.0), GeometryError);
}

TEST_CASE("offset of a cylinder is a coaxial cylinder", "[offset]") {
  const double R = 1.5, H = 3.0;
  const auto s = extrude(make_circle({0, 0, 0}, ez, R), {0, 0, H});
  const double d = 0.4;
  const double sgn = orientation(s, 0.3, 0.5, {0, 0, s.evaluate(0.3, 0.5).z});
  const auto o = offset_surface(s, d);
  CHECK(surface_deviation(s, o, d) <= 1e-7);
  CHECK(radial_error(o, [&](const Vec3& p) { return std::hypot(p.x, p.y) - (R + sgn * d); }) <= 1e-7);
  CHECK(o.num_v() == 6);  // straight along the axis: a single span
  CHECK_THROWS_AS(offset_surface(s, -sgn * 1.6), GeometryError);
}

TEST_CASE("offset of a torus is a torus with tube radius r + d", "[offset]") {
  const double Rm = 3.0, r = 1.0, d = 0.3, tol = 1e-6;
  const auto s = revolve(make_circle({Rm, 0, 0}, {0, 1, 0}, r), {0, 0, 0}, ez, 2 * pi);
  // distance of p from the tube centre circle (radius Rm in the xy-plane)
  const auto tube = [&](const Vec3& p) { return std::hypot(std::hypot(p.x, p.y) - Rm, p.z); };
  const Vec3 p0 = s.evaluate(0.3, 0.3), c0 = Rm * normalized({p0.x, p0.y, 0});
  const double sgn = orientation(s, 0.3, 0.3, c0);
  const auto o = offset_surface(s, d, tol);
  CHECK(surface_deviation(s, o, d) <= tol);
  CHECK(radial_error(o, [&](const Vec3& p) { return tube(p) - (r + sgn * d); }) <= tol);
  CHECK_THROWS_AS(offset_surface(s, -sgn * 1.2), GeometryError);  // beyond the tube radius
  // offsetting away from the tube centre is limited by the inner equator (radius of curvature
  // (Rm - r) / 1 in the direction around the axis, i.e. Rm - r = 2 on the inside)
  CHECK_THROWS_AS(offset_surface(s, sgn * 2.5), GeometryError);
}

TEST_CASE("offset of a plane is a parallel plane", "[offset]") {
  const auto s = make_rectangle({1, 2, 3}, {2, 0, 0}, {0, 1, 1});
  const Vec3 n = normalized(cross(Vec3{2, 0, 0}, Vec3{0, 1, 1}));
  const auto o = offset_surface(s, -0.75);
  CHECK(o.num_u() == 6);
  CHECK(o.num_v() == 6);
  for (int i = 0; i <= 10; ++i)
    for (int j = 0; j <= 10; ++j) {
      const double u = i / 10.0, v = j / 10.0;
      CHECK(distance(o.evaluate(u, v), s.evaluate(u, v) - 0.75 * n) <= 1e-12);
    }
  CHECK_NOTHROW(offset_surface(s, 1e6));  // a plane has no curvature limit
  CHECK_THROWS_AS(offset_surface(s, 0.1, 0.0), GeometryError);
}

TEST_CASE("offset of a freeform patch matches the brute-force distance", "[offset]") {
  std::vector<Vec3> P;
  for (int i = 0; i < 5; ++i)
    for (int j = 0; j < 4; ++j)
      P.push_back({1.0 * i, 1.0 * j, 0.3 * std::sin(1.3 * i + 0.7 * j)});
  const NurbsSurface s{3, 3, clamped_uniform_knots(3, 5), clamped_uniform_knots(3, 4), 5, 4, P};
  const double d = 0.2, tol = 1e-6;
  const auto o = offset_surface(s, d, tol);
  CHECK(surface_deviation(s, o, d) <= tol);
  for (int i = 1; i < 10; ++i)
    for (int j = 1; j < 10; ++j) {
      const Vec3 p = o.evaluate(i / 10.0, j / 10.0);
      CHECK(std::abs(closest_point(s, p).distance - d) <= 2 * tol);
    }
}

