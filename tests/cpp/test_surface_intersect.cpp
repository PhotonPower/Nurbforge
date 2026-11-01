#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <functional>
#include <numbers>
#include <vector>

#include "nurbforge/analysis/properties.hpp"
#include "nurbforge/analysis/surface_intersect.hpp"
#include "nurbforge/geometry/primitives.hpp"
#include "nurbforge/geometry/transform.hpp"
#include "nurbforge/modeling/construct.hpp"

using namespace nf;
constexpr double pi = std::numbers::pi;

namespace {

NurbsSurface sphere(const Vec3& c, double R) {
  const auto meridian = make_arc(c, {0, 0, -1}, {1, 0, 0}, R, 0.0, pi);
  return revolve(meridian, c, {0, 0, 1}, 2 * pi);
}

NurbsSurface cylinder(const Vec3& base, const Vec3& axis, double R, double len) {
  return extrude(make_circle(base, axis, R), len * normalized(axis));
}

double total_length(const SurfaceIntersection& r) {
  double L = 0.0;
  for (const auto& c : r.curves) L += arc_length(c.curve);
  return L;
}

// Every traced point and the fitted curves satisfy both surface equations; the p-curves map
// onto the 3D curve.
void check_curves(const SurfaceIntersection& r, const NurbsSurface& a, const NurbsSurface& b,
                  const std::function<double(const Vec3&)>& fa, const std::function<double(const Vec3&)>& fb,
                  double fit_tol) {
  for (const auto& c : r.curves) {
    for (const auto& p : c.points) {
      CHECK(std::abs(fa(p.point)) < 1e-9);
      CHECK(std::abs(fb(p.point)) < 1e-9);
      CHECK(distance(a.evaluate(p.u1, p.v1), p.point) < 1e-9);
      CHECK(distance(b.evaluate(p.u2, p.v2), p.point) < 1e-9);
    }
    const double t0 = c.curve.domain_start(), t1 = c.curve.domain_end();
    for (int i = 0; i <= 50; ++i) {
      const double t = t0 + (t1 - t0) * i / 50.0;
      const Vec3 q = c.curve.evaluate(t);
      CHECK(std::abs(fa(q)) < fit_tol);
      CHECK(std::abs(fb(q)) < fit_tol);
      const Vec3 w1 = c.pcurve1.evaluate(t), w2 = c.pcurve2.evaluate(t);
      CHECK(distance(a.evaluate(w1.x, w1.y), q) < fit_tol);
      CHECK(distance(b.evaluate(w2.x, w2.y), q) < fit_tol);
    }
  }
}

}  // namespace

TEST_CASE("plane / sphere: a circle", "[ssi]") {
  const double R = 2.0, h = 0.8;
  const auto s = sphere({0, 0, 0}, R);
  const auto plane = make_rectangle({-3, -3, h}, {6, 0, 0}, {0, 6, 0});
  const auto r = intersect(s, plane);
  REQUIRE(!r.curves.empty());
  CHECK(r.points.empty());
  const double rho = std::sqrt(R * R - h * h);
  CHECK(std::abs(total_length(r) - 2 * pi * rho) < 1e-5);
  check_curves(r, s, plane, [&](const Vec3& p) { return norm(p) - R; },
               [&](const Vec3& p) { return p.z - h; }, 1e-6);
}

TEST_CASE("sphere / sphere: a circle", "[ssi]") {
  const auto a = sphere({0, 0, 0}, 2.0), b = sphere({1.5, 0.3, 0.4}, 1.6);
  const auto r = intersect(a, b);
  REQUIRE(!r.curves.empty());
  // intersection circle of two spheres: radius from the distance of the centres
  const double d = norm(Vec3{1.5, 0.3, 0.4});
  const double x = (d * d + 4.0 - 1.6 * 1.6) / (2 * d);
  CHECK(std::abs(total_length(r) - 2 * pi * std::sqrt(4.0 - x * x)) < 1e-5);
  check_curves(r, a, b, [](const Vec3& p) { return norm(p) - 2.0; },
               [](const Vec3& p) { return distance(p, {1.5, 0.3, 0.4}) - 1.6; }, 1e-6);
}

TEST_CASE("plane through a torus: two circles", "[ssi]") {
  const double R = 5.0, rr = 1.5;
  const auto torus = revolve(make_circle({R, 0, 0}, {0, 1, 0}, rr), {0, 0, 0}, {0, 0, 1}, 2 * pi);
  const auto plane = make_rectangle({-8, -8, 0}, {16, 0, 0}, {0, 16, 0});
  const auto r = intersect(torus, plane);
  CHECK(std::abs(total_length(r) - 2 * pi * (R + rr) - 2 * pi * (R - rr)) < 1e-5);
  check_curves(r, torus, plane,
               [&](const Vec3& p) { return std::hypot(std::hypot(p.x, p.y) - R, p.z) - rr; },
               [](const Vec3& p) { return p.z; }, 1e-6);
}

TEST_CASE("two cylinders of different radii: two closed space curves", "[ssi]") {
  const auto big = cylinder({0, 0, -4}, {0, 0, 1}, 2.0, 8.0);
  const auto small = cylinder({-4, 0, 0.3}, {1, 0, 0}, 1.0, 8.0);
  const auto r = intersect(big, small);
  REQUIRE(r.curves.size() >= 2);
  check_curves(r, big, small, [](const Vec3& p) { return std::hypot(p.x, p.y) - 2.0; },
               [](const Vec3& p) { return std::hypot(p.y, p.z - 0.3) - 1.0; }, 1e-6);
  // every point has |x| > 0 on one of the two sides: both loops are present
  bool left = false, right = false;
  for (const auto& c : r.curves)
    for (const auto& p : c.points) (p.point.x < 0 ? left : right) = true;
  CHECK(left);
  CHECK(right);
}

TEST_CASE("plane / wavy patch: open curves ending on the boundary", "[ssi]") {
  std::vector<Vec3> P;
  for (int i = 0; i < 6; ++i)
    for (int j = 0; j < 5; ++j) P.push_back({double(i), double(j), std::sin(1.7 * i) * std::cos(1.3 * j)});
  const NurbsSurface patch{3, 3, clamped_uniform_knots(3, 6), clamped_uniform_knots(3, 5), 6, 5, P};
  const auto plane = make_rectangle({-1, -1, 0.1}, {7, 0, 0}, {0, 6, 0});
  const auto r = intersect(patch, plane);
  REQUIRE(!r.curves.empty());
  check_curves(r, patch, plane, [&](const Vec3&) { return 0.0; }, [](const Vec3& p) { return p.z - 0.1; }, 1e-6);
  for (const auto& c : r.curves) {
    if (c.closed) continue;
    for (const auto* e : {&c.points.front(), &c.points.back()}) {
      const bool on_edge = e->u1 <= 1e-9 || e->u1 >= 1 - 1e-9 || e->v1 <= 1e-9 || e->v1 >= 1 - 1e-9;
      CHECK(on_edge);
    }
  }
}

TEST_CASE("touching and disjoint surfaces", "[ssi]") {
  const auto s = sphere({0, 0, 0}, 1.0);
  const auto touching = make_rectangle({-2, -2, 1}, {4, 0, 0}, {0, 4, 0});
  const auto r = intersect(s, touching);
  CHECK(r.curves.empty());
  REQUIRE(r.points.size() == 1);
  CHECK(distance(r.points[0].point, {0, 0, 1}) < 1e-5);
  const auto far = make_rectangle({-2, -2, 1.5}, {4, 0, 0}, {0, 4, 0});
  const auto none = intersect(s, far);
  CHECK(none.curves.empty());
  CHECK(none.points.empty());
}

TEST_CASE("surface intersection: tangential contacts are points", "[ssi][tangent]") {
  const auto sph = [](const Vec3& c, double r) {
    return revolve(make_arc(c, {0, 0, -1}, {1, 0, 0}, r, 0.0, pi), c, {0, 0, 1}, 2 * pi);
  };
  // internal contact on the seam of both spheres
  const auto r1 = intersect(sph({0, 0, 0}, 1.0), sph({0.6, 0, 0}, 0.4));
  CHECK(r1.curves.empty());
  REQUIRE_FALSE(r1.points.empty());
  CHECK(distance(r1.points.front().point, Vec3{1, 0, 0}) < 1e-5);
  // two cylinders touching along a line
  const auto c1 = extrude(make_circle({0, 0, 0}, {0, 0, 1}, 0.7), {0, 0, 2});
  const auto c2 = extrude(make_circle({1.2 * std::cos(0.4), 1.2 * std::sin(0.4), 0}, {0, 0, 1}, 0.5), {0, 0, 2});
  const auto r2 = intersect(c1, c2);
  CHECK(r2.curves.empty());
  CHECK_FALSE(r2.points.empty());
}

TEST_CASE("surface intersection: a curve touching the domain boundary is not cut there", "[ssi][tangent]") {
  // circle r 0.7 on the plane [-0.7, 1]^2 touches the domain edges at (-0.7, 0) and (0, -0.7)
  const auto side = extrude(make_circle({0, 0, 0}, {0, 0, 1}, 0.7), {0, 0, 2});
  const auto plane = make_rectangle({-0.7, -0.7, 0.5}, {1.7, 0, 0}, {0, 1.7, 0});
  const auto r = intersect(plane, side);
  double length = 0.0;
  for (const auto& c : r.curves) {
    for (std::size_t k = 1; k < c.points.size(); ++k) length += distance(c.points[k].point, c.points[k - 1].point);
    for (const Vec3& end : {c.points.front().point, c.points.back().point})
      for (const Vec3& touch : {Vec3{-0.7, 0, 0.5}, Vec3{0, -0.7, 0.5}}) CHECK(distance(end, touch) > 1e-3);
  }
  CHECK(std::abs(length - 2 * pi * 0.7) < 1e-3);  // polyline: chords shorter than the arc
}
