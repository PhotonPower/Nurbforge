#include <catch2/catch_test_macros.hpp>
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

// Max distance between two curves sampled on the domain of `a` (both share the domain).
double max_deviation(const NurbsCurve& a, const NurbsCurve& b, std::size_t n = 401) {
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

// Rational cubic with an interior double knot.
NurbsCurve rational_cubic() {
  return {3,
          {0, 0, 0, 0, 0.3, 0.3, 0.7, 1, 1, 1, 1},
          {{0, 0, 0}, {1, 2, 0}, {2, 3, 1}, {3, 1, 1}, {4, 0, 2}, {5, 2, 2}, {6, 1, 0}},
          {1.0, 0.5, 2.0, 1.0, 0.7, 1.3, 1.0}};
}

// Unclamped uniform quadratic: domain [2, 5] of knots 0..8.
NurbsCurve unclamped_quadratic() {
  return {2,
          {0, 1, 2, 3, 4, 5, 6, 7},
          {{0, 0, 0}, {1, 1, 0}, {2, -1, 0}, {3, 1, 1}, {4, 0, 0}},
          {}};
}

}  // namespace

TEST_CASE("knot refinement keeps the shape", "[refine]") {
  const auto c = rational_cubic();
  const std::vector<double> X{0.9, 0.1, 0.3, 0.5, 0.5, 0.05, 0.7, 0.7};
  const auto r = c.refine_knots(X);
  CHECK(r.num_control_points() == c.num_control_points() + X.size());
  CHECK(r.knots().size() == c.knots().size() + X.size());
  CHECK(knot_multiplicity(r.knots(), 0.3) == 3);
  CHECK(knot_multiplicity(r.knots(), 0.5) == 2);
  CHECK(max_deviation(c, r) < kExact);
}

TEST_CASE("knot refinement equals repeated knot insertion", "[refine]") {
  const auto c = rational_cubic();
  const auto r = c.refine_knots({0.5, 0.2});
  const auto s = c.insert_knot(0.2).insert_knot(0.5);
  REQUIRE(r.num_control_points() == s.num_control_points());
  for (std::size_t i = 0; i < r.num_control_points(); ++i) {
    CHECK(distance(r.control_points()[i], s.control_points()[i]) < kExact);
    CHECK(std::abs(r.weights()[i] - s.weights()[i]) < kExact);
  }
}

TEST_CASE("knot refinement of an unclamped curve, including domain ends", "[refine]") {
  const auto c = unclamped_quadratic();
  const auto r = c.refine_knots({2.0, 3.5, 3.5, 5.0});  // ends reach multiplicity p = 2
  CHECK(knot_multiplicity(r.knots(), 2.0) == 2);
  CHECK(knot_multiplicity(r.knots(), 5.0) == 2);
  CHECK(r.domain_start() == 2.0);
  CHECK(r.domain_end() == 5.0);
  CHECK(max_deviation(c, r) < kExact);
  CHECK_THROWS_AS(c.refine_knots({2.0, 2.0}), GeometryError);
}

TEST_CASE("knot refinement of an exact circle", "[refine]") {
  const auto c = make_circle({1, 2, 3}, {0, 0, 1}, 2.0);
  const auto r = c.refine_knots({0.1, 0.33, 0.6, 0.6, 0.95});
  CHECK(max_deviation(c, r) < kExact);
  for (double u = 0.0; u <= 1.0; u += 0.01) CHECK(std::abs(distance(r.evaluate(u), {1, 2, 3}) - 2.0) < kExact);
}

TEST_CASE("knot refinement rejects invalid input", "[refine]") {
  const auto c = rational_cubic();
  CHECK_THROWS_AS(c.refine_knots({1.5}), GeometryError);
  CHECK_THROWS_AS(c.refine_knots({-0.1}), GeometryError);
  CHECK_THROWS_AS(c.refine_knots({0.3, 0.3}), GeometryError);  // 0.3 already double
  CHECK_THROWS_AS(c.refine_knots({0.0}), GeometryError);       // clamped end has p+1
  CHECK_THROWS_AS(c.refine_knots({0.5, 0.5, 0.5, 0.5}), GeometryError);
  CHECK(c.refine_knots({}).num_control_points() == c.num_control_points());
}

TEST_CASE("Bezier decomposition of a rational cubic", "[refine][bezier]") {
  const auto c = rational_cubic();
  const auto segs = c.decompose_bezier();
  REQUIRE(segs.size() == 3);  // spans [0,.3], [.3,.7], [.7,1]
  const double bounds[] = {0.0, 0.3, 0.7, 1.0};
  for (std::size_t k = 0; k < segs.size(); ++k) {
    const auto& s = segs[k];
    CHECK(s.degree() == 3);
    CHECK(s.num_control_points() == 4);
    CHECK(s.domain_start() == bounds[k]);
    CHECK(s.domain_end() == bounds[k + 1]);
    CHECK(max_deviation(s, c) < kExact);  // sampled on the segment's own interval
  }
  // C0 joints: consecutive segments share the end control point
  for (std::size_t k = 0; k + 1 < segs.size(); ++k)
    CHECK(distance(segs[k].control_points().back(), segs[k + 1].control_points().front()) < kExact);
  CHECK(distance(segs.front().control_points().front(), c.control_points().front()) < kExact);
  CHECK(distance(segs.back().control_points().back(), c.control_points().back()) < kExact);
}

TEST_CASE("Bezier decomposition of an unclamped curve", "[refine][bezier]") {
  const auto c = unclamped_quadratic();
  const auto segs = c.decompose_bezier();
  REQUIRE(segs.size() == 3);  // spans [2,3], [3,4], [4,5]
  for (std::size_t k = 0; k < segs.size(); ++k) {
    CHECK(segs[k].domain_start() == 2.0 + static_cast<double>(k));
    CHECK(max_deviation(segs[k], c) < kExact);
  }
}

TEST_CASE("Bezier decomposition of a circle gives four quarter arcs", "[refine][bezier]") {
  const auto c = make_circle({0, 0, 0}, {0, 0, 1}, 3.0);
  const auto segs = c.decompose_bezier();
  REQUIRE(segs.size() == 4);
  for (const auto& s : segs) {
    CHECK(s.is_rational());
    CHECK(std::abs(s.weights()[1] / s.weights()[0] - std::cos(std::numbers::pi / 4)) < kExact);
    CHECK(max_deviation(s, c) < kExact);
  }
}

TEST_CASE("Bezier segment of a single-span curve is the curve itself", "[refine][bezier]") {
  const auto c = make_line({0, 0, 0}, {1, 2, 3});
  const auto segs = c.decompose_bezier();
  REQUIRE(segs.size() == 1);
  CHECK(max_deviation(segs[0], c) < kExact);
}

TEST_CASE("surface knot refinement in u and v keeps the shape", "[refine][surface]") {
  // torus: rational in both directions
  const auto profile = make_circle({5, 0, 0}, {0, 1, 0}, 1.5);
  const auto torus = revolve(profile, {0, 0, 0}, {0, 0, 1}, 2 * std::numbers::pi);
  const auto ru = torus.refine_knots_u({0.1, 0.4, 0.4, 0.8});
  CHECK(ru.num_u() == torus.num_u() + 4);
  CHECK(ru.num_v() == torus.num_v());
  CHECK(max_deviation(torus, ru) < kExact);
  const auto rv = torus.refine_knots_v({0.2, 0.6});
  CHECK(rv.num_v() == torus.num_v() + 2);
  CHECK(max_deviation(torus, rv) < kExact);
  CHECK_THROWS_AS(torus.refine_knots_v({0.5, 0.5}), GeometryError);  // 0.5 already double
}

TEST_CASE("surface Bezier decomposition", "[refine][bezier][surface]") {
  const auto profile = make_circle({5, 0, 0}, {0, 1, 0}, 1.5);
  const auto torus = revolve(profile, {0, 0, 0}, {0, 0, 1}, 3 * std::numbers::pi / 2);
  const auto patches = torus.decompose_bezier();
  REQUIRE(patches.size() == 3 * 4);  // 3 spans around the axis (270°), 4 around the tube
  for (const auto& p : patches) {
    CHECK(p.num_u() == 3);
    CHECK(p.num_v() == 3);
    CHECK(max_deviation(p, torus, 11) < kExact);
  }
  // row-major order: u span outer
  CHECK(patches[0].domain_u()[0] == 0.0);
  CHECK(patches[1].domain_u()[0] == 0.0);
  CHECK(patches[1].domain_v()[0] == patches[0].domain_v()[1]);
  CHECK(patches[4].domain_u()[0] == patches[0].domain_u()[1]);
}
