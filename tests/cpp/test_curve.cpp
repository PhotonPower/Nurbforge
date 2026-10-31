#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <numbers>

#include "nurbforge/geometry/primitives.hpp"
#include "nurbforge/nurbs/curve.hpp"

using Catch::Approx;
using namespace nf;

TEST_CASE("line evaluation", "[curve]") {
  const auto c = make_line({0, 0, 0}, {2, 4, 6});
  const auto p = c.evaluate(0.5);
  CHECK(p.x == Approx(1)); CHECK(p.y == Approx(2)); CHECK(p.z == Approx(3));
}

TEST_CASE("exact circle stays on radius", "[curve]") {
  const auto c = make_circle({1, 2, 3}, {0, 0, 1}, 2.5);
  CHECK(c.num_control_points() == 9);
  CHECK(c.is_rational());
  for (double u = 0; u <= 1.0; u += 0.01) {
    const auto p = c.evaluate(u);
    CHECK(distance(p, {1, 2, 3}) == Approx(2.5).epsilon(1e-12));
    CHECK(p.z == Approx(3.0));
  }
}

TEST_CASE("uniformly scaled weights are not rational", "[curve][tolerance]") {
  const KnotVector U{0, 0, 0, 1, 1, 1};
  const std::vector<Vec3> P{{0, 0, 0}, {1, 1, 0}, {2, 0, 0}};
  CHECK_FALSE(NurbsCurve(2, U, P, {1e6, 1e6, 1e6}).is_rational());
  CHECK_FALSE(NurbsCurve(2, U, P, {1e-12, 1e-12, 1e-12}).is_rational());
  // differences below the absolute tolerance must still count when the weights are tiny
  CHECK(NurbsCurve(2, U, P, {1e-12, 2e-12, 1e-12}).is_rational());
}

TEST_CASE("arc end points", "[curve]") {
  const auto a = make_arc({0, 0, 0}, {1, 0, 0}, {0, 1, 0}, 1.0, 0.0, std::numbers::pi / 3);
  const auto e = a.evaluate(a.domain_end());
  CHECK(e.x == Approx(0.5)); CHECK(e.y == Approx(std::sqrt(3.0) / 2));
}

TEST_CASE("knot insertion preserves geometry", "[curve]") {
  const auto c = NurbsCurve::from_points(3, {{0, 0, 0}, {1, 2, 0}, {3, 2, 1}, {4, 0, 0}, {6, 1, 2}});
  const auto ci = c.insert_knot(0.3, 2);
  CHECK(ci.num_control_points() == c.num_control_points() + 2);
  for (double u = 0; u <= 1.0; u += 0.05) CHECK(distance(c.evaluate(u), ci.evaluate(u)) < 1e-12);

  const auto circ = make_circle({0, 0, 0}, {0, 0, 1}, 1.0);
  const auto circ2 = circ.insert_knot(0.1);
  for (double u = 0; u <= 1.0; u += 0.05) CHECK(distance(circ.evaluate(u), circ2.evaluate(u)) < 1e-12);
}

TEST_CASE("rational derivative of a circle is tangent", "[curve]") {
  const auto c = make_circle({0, 0, 0}, {0, 0, 1}, 1.0);
  for (double u : {0.05, 0.3, 0.61, 0.9}) {
    const auto d = c.derivatives(u, 2);
    CHECK(dot(d[0], d[1]) == Approx(0.0).margin(1e-10));
    // finite difference check
    const double h = 1e-6;
    const Vec3 fd = (c.evaluate(u + h) - c.evaluate(u - h)) / (2 * h);
    CHECK(distance(fd, d[1]) < 1e-5);
  }
}

TEST_CASE("reversal keeps the point set", "[curve]") {
  const auto c = NurbsCurve::from_points(2, {{0, 0, 0}, {1, 1, 0}, {2, 0, 0}, {3, 1, 0}});
  const auto r = c.reversed();
  CHECK(distance(c.evaluate(0.2), r.evaluate(0.8)) < 1e-12);
}
