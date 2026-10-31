#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <numbers>
#include <vector>

#include "nurbforge/analysis/properties.hpp"
#include "nurbforge/core/error.hpp"
#include "nurbforge/geometry/primitives.hpp"
#include "nurbforge/geometry/transform.hpp"
#include "nurbforge/modeling/construct.hpp"

using namespace nf;
constexpr double pi = std::numbers::pi;

namespace {

bool rel_close(double a, double b, double rel) { return std::abs(a - b) <= rel * std::max(std::abs(b), 1e-300); }

// (t, t^2, t^3) on [0, 1] as a cubic Bezier curve.
NurbsCurve twisted_cubic() {
  return {3, {0, 0, 0, 0, 1, 1, 1, 1}, {{0, 0, 0}, {1.0 / 3, 0, 0}, {2.0 / 3, 1.0 / 3, 0}, {1, 1, 1}}, {}};
}

NurbsSurface sphere(double R) {
  const auto meridian = make_arc({0, 0, 0}, {0, 0, -1}, {1, 0, 0}, R, 0.0, pi);
  return revolve(meridian, {0, 0, 0}, {0, 0, 1}, 2 * pi);
}

}  // namespace

// ------------------------------------------------------------------ curves

TEST_CASE("curvature and torsion of the twisted cubic", "[properties]") {
  const auto c = twisted_cubic();
  for (double t = 0.0; t <= 1.0; t += 0.05) {
    const double q = 36 * t * t * t * t + 36 * t * t + 4;
    const double speed2 = 1 + 4 * t * t + 9 * t * t * t * t;
    CHECK(rel_close(curvature(c, t), std::sqrt(q) / std::pow(speed2, 1.5), 1e-12));
    CHECK(rel_close(torsion(c, t), 3.0 / (9 * t * t * t * t + 9 * t * t + 1), 1e-12));
  }
}

TEST_CASE("curvature of circle and line", "[properties]") {
  const auto circle = make_circle({1, 2, 3}, {1, 1, 0}, 3.0);
  for (double u = 0.0; u <= 1.0; u += 0.03) {
    CHECK(rel_close(curvature(circle, u), 1.0 / 3.0, 1e-12));
    CHECK(std::abs(torsion(circle, u)) < 1e-12);  // planar
  }
  const auto line = make_line({0, 0, 0}, {1, 2, 3});
  CHECK(curvature(line, 0.4) == 0.0);
  CHECK_THROWS_AS(torsion(line, 0.4), GeometryError);
}

TEST_CASE("arc length", "[properties]") {
  const double R = 3.0;
  const auto circle = make_circle({0, 0, 0}, {0, 0, 1}, R);
  CHECK(rel_close(arc_length(circle), 2 * pi * R, 1e-12));
  CHECK(rel_close(arc_length(circle, 0.0, 0.25), pi * R / 2, 1e-12));  // quarter arc
  CHECK(rel_close(arc_length(circle, 0.0, 0.37) + arc_length(circle, 0.37, 1.0), 2 * pi * R, 1e-12));
  CHECK(arc_length(circle, 0.6, 0.2) == -arc_length(circle, 0.2, 0.6));
  CHECK(arc_length(circle, 0.5, 0.5) == 0.0);

  const auto pl = make_polyline({{0, 0, 0}, {3, 4, 0}, {3, 4, 12}});
  CHECK(rel_close(arc_length(pl), 5.0 + 12.0, 1e-13));

  // twisted cubic against a fine composite reference
  const auto c = twisted_cubic();
  double ref = 0.0;
  const int n = 200000;
  for (int i = 0; i < n; ++i) ref += distance(c.evaluate(double(i) / n), c.evaluate(double(i + 1) / n));
  CHECK(rel_close(arc_length(c), ref, 1e-9));
  CHECK_THROWS_AS(arc_length(c, -0.5, 0.5), GeometryError);
}

TEST_CASE("parameter at arc length", "[properties]") {
  const auto c = twisted_cubic();
  const double L = arc_length(c);
  for (double f : {0.1, 0.33, 0.5, 0.9}) {
    const double u = parameter_at_length(c, f * L);
    CHECK(rel_close(arc_length(c, 0.0, u), f * L, 1e-9));
  }
  CHECK(parameter_at_length(c, 0.0) == 0.0);
  CHECK(parameter_at_length(c, L) == 1.0);
  CHECK_THROWS_AS(parameter_at_length(c, 1.1 * L), GeometryError);
  // a circle has constant speed only per quarter (rational parametrisation): still exact
  const auto circle = make_circle({0, 0, 0}, {0, 0, 1}, 1.0);
  const double u = parameter_at_length(circle, pi / 3);
  CHECK(rel_close(arc_length(circle, 0.0, u), pi / 3, 1e-10));
}

// ------------------------------------------------------------------ surfaces

TEST_CASE("sphere curvature, also at the poles", "[properties][surface]") {
  const double R = 2.0;
  const auto s = sphere(R);
  for (double u : {0.0, 0.3, 0.75})
    for (double v : {0.0, 0.2, 0.5, 0.9, 1.0}) {
      const auto k = curvature(s, u, v);
      CHECK(rel_close(k.gaussian, 1.0 / (R * R), 1e-6));
      CHECK(rel_close(std::abs(k.mean), 1.0 / R, 1e-6));
      CHECK(std::abs(k.k1 - k.k2) < 1e-6);
      CHECK(std::abs(dot(k.dir1, k.dir2)) < 1e-9);
      CHECK(std::abs(dot(k.dir1, k.normal)) < 1e-9);
    }
}

TEST_CASE("torus curvature and principal directions", "[properties][surface]") {
  const double R = 5.0, r = 1.5;
  const auto torus = revolve(make_circle({R, 0, 0}, {0, 1, 0}, r), {0, 0, 0}, {0, 0, 1}, 2 * pi);
  for (double u : {0.1, 0.45, 0.8})
    for (double v : {0.05, 0.3, 0.55, 0.85}) {
      const Vec3 p = torus.evaluate(u, v);
      const double rho = std::hypot(p.x, p.y);
      const double cos_t = (rho - R) / r;
      const double k_ring = cos_t / rho, k_tube = 1.0 / r;
      const auto k = curvature(torus, u, v);
      CHECK(std::abs(k.gaussian - k_ring * k_tube) < 1e-10);
      CHECK(std::abs(std::abs(k.mean) - 0.5 * std::abs(k_ring + k_tube)) < 1e-10);
      // the tube curvature belongs to the meridian direction, orthogonal to the parallel
      const Vec3 parallel{-p.y / rho, p.x / rho, 0.0};
      const Vec3& tube_dir = std::abs(std::abs(k.k1) - k_tube) < 1e-8 ? k.dir1 : k.dir2;
      CHECK(std::abs(dot(tube_dir, parallel)) < 1e-8);
    }
}

TEST_CASE("cylinder and plane curvature", "[properties][surface]") {
  const auto cyl = extrude(make_circle({0, 0, 0}, {0, 0, 1}, 4.0), {0, 0, 10});
  const auto k = curvature(cyl, 0.3, 0.6);
  CHECK(std::abs(k.gaussian) < 1e-12);
  CHECK(rel_close(std::abs(k.mean), 1.0 / 8.0, 1e-12));
  const auto plane = make_rectangle({0, 0, 0}, {1, 0, 0}, {0, 2, 0});
  const auto kp = curvature(plane, 0.5, 0.5);
  CHECK(kp.gaussian == 0.0);
  CHECK(kp.k1 == 0.0);
  CHECK(std::abs(dot(kp.dir1, kp.dir2)) < 1e-15);
}

TEST_CASE("surface area", "[properties][surface]") {
  const double R = 2.0;
  const auto s = sphere(R);
  CHECK(rel_close(area(s), 4 * pi * R * R, 1e-10));
  CHECK(rel_close(area(s, 0.0, 1.0, 0.0, 0.5), 2 * pi * R * R, 1e-10));  // southern hemisphere
  const auto torus = revolve(make_circle({5, 0, 0}, {0, 1, 0}, 1.5), {0, 0, 0}, {0, 0, 1}, 2 * pi);
  CHECK(rel_close(area(torus), 4 * pi * pi * 5 * 1.5, 1e-10));
  const auto cyl = extrude(make_circle({0, 0, 0}, {0, 0, 1}, 4.0), {0, 0, 10});
  CHECK(rel_close(area(cyl), 2 * pi * 4 * 10, 1e-10));
  CHECK(rel_close(area(make_rectangle({0, 0, 0}, {2, 0, 0}, {0, 3, 0})), 6.0, 1e-14));
  // scaling by 3 multiplies the area by 9
  const auto big = torus.transformed(Transform::scaling({0, 0, 0}, 3.0));
  CHECK(rel_close(area(big), 9 * area(torus), 1e-10));
  CHECK_THROWS_AS(area(s, 0.5, 0.4, 0.0, 1.0), GeometryError);
}
