#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <numbers>

#include "nurbforge/geometry/primitives.hpp"
#include "nurbforge/modeling/construct.hpp"
#include "nurbforge/modeling/pattern.hpp"

using Catch::Approx;
using namespace nf;

TEST_CASE("revolving a line gives an exact cylinder", "[construct]") {
  const auto line = make_line({2, 0, 0}, {2, 0, 5});
  const auto cyl = revolve(line, {0, 0, 0}, {0, 0, 1}, 2 * std::numbers::pi);
  for (double u = 0; u <= 1; u += 0.05)
    for (double v = 0; v <= 1; v += 0.25) {
      const auto p = cyl.evaluate(u, v);
      CHECK(std::hypot(p.x, p.y) == Approx(2.0));
      CHECK(p.z == Approx(5 * v).margin(1e-12));
    }
}

TEST_CASE("revolving a semicircle gives a sphere", "[construct]") {
  const auto arc = make_arc({0, 0, 0}, {0, 0, -1}, {1, 0, 0}, 3.0, 0, std::numbers::pi);
  const auto sph = revolve(arc, {0, 0, 0}, {0, 0, 1}, 2 * std::numbers::pi);
  for (double u = 0; u <= 1; u += 0.1)
    for (double v = 0; v <= 1; v += 0.1) CHECK(norm(sph.evaluate(u, v)) == Approx(3.0));
}

TEST_CASE("partial revolution", "[construct]") {
  const auto line = make_line({1, 0, 0}, {1, 0, 1});
  const auto s = revolve(line, {0, 0, 0}, {0, 0, 1}, std::numbers::pi / 3);
  const auto p = s.evaluate(1.0, 0.0);
  CHECK(p.x == Approx(0.5)); CHECK(p.y == Approx(std::sqrt(3.0) / 2));
}

TEST_CASE("extrusion", "[construct]") {
  const auto c = make_circle({0, 0, 0}, {0, 0, 1}, 1.0);
  const auto s = extrude(c, {0, 0, 4});
  const auto p = s.evaluate(0.3, 0.5);
  CHECK(std::hypot(p.x, p.y) == Approx(1.0)); CHECK(p.z == Approx(2.0));
}

TEST_CASE("ruled surface between incompatible knot vectors", "[construct]") {
  const auto a = NurbsCurve::from_points(2, {{0, 0, 0}, {1, 1, 0}, {2, 0, 0}});
  const auto b = NurbsCurve::from_points(2, {{0, 0, 1}, {1, 1, 1}, {2, 2, 1}, {3, 0, 1}});
  const auto s = ruled(a, b);
  CHECK(distance(s.evaluate(0.4, 0.0), a.evaluate(0.4)) < 1e-12);
  CHECK(distance(s.evaluate(0.4, 1.0), b.evaluate(0.4)) < 1e-12);
}

TEST_CASE("polar pattern", "[pattern]") {
  const auto ts = polar_pattern({0, 0, 0}, {0, 0, 1}, 4, 2 * std::numbers::pi);
  REQUIRE(ts.size() == 4);
  const auto p = ts[1].apply({1, 0, 0});
  CHECK(p.y == Approx(1.0));
}
