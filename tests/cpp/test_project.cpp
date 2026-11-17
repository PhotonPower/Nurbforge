#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <numbers>

#include "nurbforge/core/error.hpp"
#include "nurbforge/geometry/primitives.hpp"
#include "nurbforge/modeling/construct.hpp"
#include "nurbforge/modeling/project.hpp"

using namespace nf;
constexpr double pi = std::numbers::pi;

TEST_CASE("project curves onto surfaces", "[project]") {
  SECTION("circle along -z onto a saddle z = x y") {
    const auto saddle = make_bilinear({-1, -1, 1}, {1, -1, -1}, {-1, 1, -1}, {1, 1, 1});
    const auto circle = make_circle({0, 0, 3}, {0, 0, 1}, 0.5);
    const auto pr = project_curve(circle, saddle, {0, 0, -1}, 1e-7);
    for (int i = 0; i <= 50; ++i) {
      const double t = i / 50.0;
      const Vec3 c = circle.evaluate(t), p = pr.curve.evaluate(t);
      CHECK(std::abs(p.x - c.x) <= 2e-7);
      CHECK(std::abs(p.y - c.y) <= 2e-7);
      CHECK(std::abs(p.z - p.x * p.y) <= 2e-7);  // on the saddle
      const Vec3 w = pr.pcurve.evaluate(t);
      CHECK(distance(saddle.evaluate(w.x, w.y), p) <= 2e-7);
    }
  }
  SECTION("normal projection: a line onto a cylinder lies on the cylinder") {
    const auto cyl = extrude(make_circle({0, 0, 0}, {0, 0, 1}, 1.0), {0, 0, 4});
    const auto line = make_line({0.0, -2.0, 0.5}, {0.0, -2.5, 3.5});  // outside, at angle -90
    const auto pr = project_curve(line, cyl, {0, 0, 0}, 1e-7);
    for (int i = 0; i <= 20; ++i) {
      const Vec3 p = pr.curve.evaluate(i / 20.0);
      CHECK(std::abs(std::hypot(p.x, p.y) - 1.0) <= 2e-7);
      CHECK(std::abs(p.x) <= 2e-7);  // the generator at angle -90 degrees
    }
  }
  SECTION("errors: ray missing the surface") {
    const auto plate = make_rectangle({0, 0, 0}, {1, 0, 0}, {0, 1, 0});
    const auto far = make_line({5, 5, 1}, {6, 5, 1});
    CHECK_THROWS_AS(project_curve(far, plate, {0, 0, -1}), GeometryError);
  }
}
