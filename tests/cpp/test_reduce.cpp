#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <algorithm>
#include <cmath>
#include <numbers>
#include <vector>

#include "nurbforge/core/error.hpp"
#include "nurbforge/geometry/primitives.hpp"
#include "nurbforge/modeling/construct.hpp"
#include "nurbforge/nurbs/curve.hpp"
#include "nurbforge/nurbs/surface.hpp"

using namespace nf;

namespace {

constexpr double kExact = 1e-12;

double max_deviation(const NurbsCurve& a, const NurbsCurve& b, std::size_t n = 801) {
  double d = 0.0;
  const double u0 = a.domain_start(), u1 = a.domain_end();
  for (std::size_t i = 0; i < n; ++i) {
    const double u = u0 + (u1 - u0) * static_cast<double>(i) / static_cast<double>(n - 1);
    d = std::max(d, distance(a.evaluate(u), b.evaluate(u)));
  }
  return d;
}

double max_deviation(const NurbsSurface& a, const NurbsSurface& b, std::size_t n = 41) {
  double d = 0.0;
  const auto du = a.domain_u(), dv = a.domain_v();
  for (std::size_t i = 0; i < n; ++i)
    for (std::size_t j = 0; j < n; ++j) {
      const double u = du[0] + (du[1] - du[0]) * static_cast<double>(i) / static_cast<double>(n - 1);
      const double v = dv[0] + (dv[1] - dv[0]) * static_cast<double>(j) / static_cast<double>(n - 1);
      d = std::max(d, distance(a.evaluate(u, v), b.evaluate(u, v)));
    }
  return d;
}

void check_same_polygon(const NurbsCurve& a, const NurbsCurve& b, double eps = 1e-10) {
  REQUIRE(a.degree() == b.degree());
  REQUIRE(a.knots().size() == b.knots().size());
  REQUIRE(a.num_control_points() == b.num_control_points());
  for (std::size_t i = 0; i < a.knots().size(); ++i) CHECK(std::abs(a.knots()[i] - b.knots()[i]) < eps);
  // compare homogeneous points: weights may differ by a common factor only if both are scaled,
  // which the algorithms here never do
  for (std::size_t i = 0; i < a.num_control_points(); ++i) {
    CHECK(distance(a.control_points()[i], b.control_points()[i]) < eps);
    CHECK(std::abs(a.weights()[i] - b.weights()[i]) < eps);
  }
}

NurbsCurve rational_cubic() {
  return {3,
          {0, 0, 0, 0, 0.3, 0.3, 0.7, 1, 1, 1, 1},
          {{0, 0, 0}, {1, 2, 0}, {2, 3, 1}, {3, 1, 1}, {4, 0, 2}, {5, 2, 2}, {6, 1, 0}},
          {1.0, 0.5, 2.0, 1.0, 0.7, 1.3, 1.0}};
}

// A genuine cubic B-spline (not an elevated quadratic) with simple interior knots.
NurbsCurve wavy_cubic() {
  return NurbsCurve::from_points(
      3, {{0, 0, 0}, {1, 3, 0}, {2, -2, 0}, {3, 3, 1}, {4, -1, 0}, {5, 2, -1}, {6, 0, 0}});
}

}  // namespace

// ------------------------------------------------------------------ knot removal

TEST_CASE("inserted knots are removed exactly", "[remove]") {
  const auto c = rational_cubic();
  check_same_polygon(c.insert_knot(0.5, 3).remove_knot(0.5, 3), c);
  check_same_polygon(c.insert_knot(0.3).remove_knot(0.3), c);  // back to multiplicity 2
  check_same_polygon(c.refine_knots({0.1, 0.8, 0.8}).remove_knot(0.8, 2).remove_knot(0.1), c);
}

TEST_CASE("knot removal on an exact circle", "[remove]") {
  const auto c = make_circle({0, 0, 0}, {0, 0, 1}, 5.0);
  check_same_polygon(c.refine_knots({0.1, 0.6}).remove_knot(0.6).remove_knot(0.1), c);
  // the quarter-arc joints are genuine: not removable at tight tolerance
  const auto [same, removed] = c.try_remove_knot(0.25, 1);
  CHECK(removed == 0);
  CHECK(same.num_control_points() == c.num_control_points());
  CHECK_THROWS_AS(c.remove_knot(0.25), GeometryError);
}

TEST_CASE("try_remove_knot stops at the tolerance and respects it", "[remove]") {
  const auto c = wavy_cubic();
  const double u = c.knots()[5];  // second interior knot
  CHECK(c.try_remove_knot(u, 1, 1e-10).second == 0);
  for (double tol : {0.05, 0.3, 1.0, 3.0}) {
    const auto [r, removed] = c.try_remove_knot(u, 1, tol);
    if (removed == 1) {
      CHECK(r.num_control_points() == c.num_control_points() - 1);
      CHECK(max_deviation(c, r) <= tol);
    }
  }
  CHECK(c.try_remove_knot(u, 1, 10.0).second == 1);
}

TEST_CASE("knot removal is limited by the multiplicity", "[remove]") {
  const auto c = rational_cubic().insert_knot(0.5, 2);
  const auto [r, removed] = c.try_remove_knot(0.5, 5);
  CHECK(removed == 2);
  CHECK(knot_multiplicity(r.knots(), 0.5) == 0);
  CHECK_THROWS_AS(c.remove_knot(0.5, 3), GeometryError);
}

TEST_CASE("knot removal rejects non-interior knots", "[remove]") {
  const auto c = rational_cubic();
  CHECK_THROWS_AS(c.remove_knot(0.0), GeometryError);
  CHECK_THROWS_AS(c.remove_knot(1.0), GeometryError);
  CHECK_THROWS_AS(c.remove_knot(0.42), GeometryError);
}

// ------------------------------------------------------------------ degree reduction

TEST_CASE("degree reduction undoes degree elevation exactly", "[reduce]") {
  const auto c = rational_cubic();
  check_same_polygon(c.elevate_degree(1).reduce_degree(), c);
  check_same_polygon(c.elevate_degree(2).reduce_degree().reduce_degree(), c);
  const auto circle = make_circle({1, 2, 3}, {0, 1, 1}, 2.0);
  check_same_polygon(circle.elevate_degree(1).reduce_degree(), circle);
  const auto pl = make_polyline({{0, 0, 0}, {1, 0, 0}, {1, 1, 0}});
  check_same_polygon(pl.elevate_degree(1).reduce_degree(), pl);
}

TEST_CASE("degree reduction of an unclamped elevated curve", "[reduce]") {
  const NurbsCurve c{2, {0, 1, 2, 3, 4, 5, 6, 7}, {{0, 0, 0}, {1, 1, 0}, {2, -1, 0}, {3, 1, 1}, {4, 0, 0}}, {}};
  const auto r = c.elevate_degree(1).reduce_degree();
  CHECK(r.degree() == 2);
  CHECK(max_deviation(c, r) < kExact);
  check_same_polygon(r, c.clamped());
}

TEST_CASE("degree reduction fails for a genuine cubic at tight tolerance", "[reduce]") {
  const auto c = wavy_cubic();
  CHECK_THROWS_AS(c.reduce_degree(), GeometryError);
  CHECK_THROWS_AS(make_line({0, 0, 0}, {1, 1, 1}).reduce_degree(), GeometryError);  // degree 1
}

TEST_CASE("approximate degree reduction stays within the tolerance", "[reduce]") {
  const auto c = wavy_cubic();
  CHECK_THROWS_AS(c.reduce_degree(0.5), GeometryError);
  const double tol = GENERATE(2.0, 5.0);
  const auto r = c.reduce_degree(tol);
  CHECK(r.degree() == 2);
  CHECK(max_deviation(c, r, 4001) <= tol);
  // simple knots are kept (multiplicity 1)
  for (std::size_t k = 4; k < 7; ++k) CHECK(knot_multiplicity(r.knots(), c.knots()[k]) == 1);
}

TEST_CASE("nearly reducible curve: perturbation below the tolerance", "[reduce]") {
  auto pts = rational_cubic().elevate_degree(1).control_points();
  pts[3] = pts[3] + Vec3{1e-6, 0, 0};
  const auto e = rational_cubic().elevate_degree(1);
  const NurbsCurve perturbed{e.degree(), e.knots(), pts, e.weights()};
  CHECK_THROWS_AS(perturbed.reduce_degree(1e-9), GeometryError);
  const auto r = perturbed.reduce_degree(1e-4);
  CHECK(r.degree() == 3);
  CHECK(max_deviation(perturbed, r) <= 1e-4);
}

// ------------------------------------------------------------------ surfaces

TEST_CASE("surface knot removal and degree reduction", "[reduce][surface]") {
  const auto profile = make_circle({5, 0, 0}, {0, 1, 0}, 1.5);
  const auto torus = revolve(profile, {0, 0, 0}, {0, 0, 1}, 3 * std::numbers::pi / 2);

  const auto ru = torus.refine_knots_u({0.2, 0.9}).remove_knot_u(0.2).remove_knot_u(0.9);
  CHECK(ru.num_u() == torus.num_u());
  CHECK(max_deviation(torus, ru) < kExact);
  const auto rv = torus.refine_knots_v({0.4}).remove_knot_v(0.4);
  CHECK(rv.num_v() == torus.num_v());
  CHECK_THROWS_AS(torus.remove_knot_v(0.25), GeometryError);  // genuine joint of the tube circle

  const auto du = torus.elevate_degree_u(1).reduce_degree_u();
  CHECK(du.degree_u() == 2);
  CHECK(du.num_u() == torus.num_u());
  CHECK(max_deviation(torus, du) < kExact);
  const auto dv = torus.elevate_degree_v(2).reduce_degree_v().reduce_degree_v();
  CHECK(dv.degree_v() == 2);
  CHECK(dv.num_v() == torus.num_v());
  CHECK(max_deviation(torus, dv) < kExact);
  CHECK_THROWS_AS(torus.reduce_degree_u(), GeometryError);
}
