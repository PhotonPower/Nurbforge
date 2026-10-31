#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <algorithm>
#include <cmath>
#include <numbers>
#include <vector>

#include "nurbforge/geometry/primitives.hpp"
#include "nurbforge/modeling/construct.hpp"
#include "nurbforge/nurbs/curve.hpp"
#include "nurbforge/nurbs/surface.hpp"

using namespace nf;

namespace {

constexpr double kExact = 1e-12;

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

// Rational cubic: interior knots 0.3 (multiplicity 2) and 0.7 (1) -> 3 Bezier segments.
NurbsCurve rational_cubic() {
  return {3,
          {0, 0, 0, 0, 0.3, 0.3, 0.7, 1, 1, 1, 1},
          {{0, 0, 0}, {1, 2, 0}, {2, 3, 1}, {3, 1, 1}, {4, 0, 2}, {5, 2, 2}, {6, 1, 0}},
          {1.0, 0.5, 2.0, 1.0, 0.7, 1.3, 1.0}};
}

NurbsCurve unclamped_quadratic() {
  return {2, {0, 1, 2, 3, 4, 5, 6, 7}, {{0, 0, 0}, {1, 1, 0}, {2, -1, 0}, {3, 1, 1}, {4, 0, 0}}, {}};
}

}  // namespace

TEST_CASE("degree elevation of a rational cubic keeps the shape", "[elevate]") {
  const std::size_t t = GENERATE(1, 2, 3);
  const auto c = rational_cubic();
  const auto e = c.elevate_degree(t);
  CHECK(e.degree() == 3 + t);
  CHECK(e.num_control_points() == 7 + 3 * t);  // n+1 + t * (number of segments)
  CHECK(knot_multiplicity(e.knots(), 0.3) == 2 + t);
  CHECK(knot_multiplicity(e.knots(), 0.7) == 1 + t);
  CHECK(knot_multiplicity(e.knots(), 0.0) == 3 + t + 1);
  CHECK(e.is_clamped());
  CHECK(max_deviation(c, e) < kExact);
}

TEST_CASE("degree elevation of a line gives the expected Bezier points", "[elevate]") {
  const auto l = make_line({0, 0, 0}, {3, 6, 9});
  const auto q = l.elevate_degree(2);
  REQUIRE(q.num_control_points() == 4);
  for (std::size_t i = 0; i < 4; ++i) {
    const double s = static_cast<double>(i) / 3.0;
    CHECK(distance(q.control_points()[i], Vec3{3 * s, 6 * s, 9 * s}) < kExact);
  }
}

TEST_CASE("degree elevation of a polyline (degree 1, C0 knots)", "[elevate]") {
  const auto pl = make_polyline({{0, 0, 0}, {1, 0, 0}, {1, 2, 0}, {3, 2, 1}});
  const auto e = pl.elevate_degree(2);
  CHECK(e.degree() == 3);
  CHECK(max_deviation(pl, e) < kExact);
  // corners stay corners: interior knots now have multiplicity 3 = degree
  const auto [vals, mults] = unique_knots(e.knots());
  for (std::size_t i = 1; i + 1 < vals.size(); ++i) CHECK(mults[i] == 3);
}

TEST_CASE("degree elevation keeps an exact circle exact", "[elevate]") {
  const auto c = make_circle({1, -2, 0.5}, {1, 1, 1}, 4.0);
  const auto e = c.elevate_degree(1);
  CHECK(e.degree() == 3);
  CHECK(max_deviation(c, e) < kExact);
  for (const auto& p : e.sample(200)) CHECK(std::abs(distance(p, {1, -2, 0.5}) - 4.0) < kExact);
}

TEST_CASE("clamping and elevating an unclamped curve", "[elevate]") {
  const auto c = unclamped_quadratic();
  CHECK_FALSE(c.is_clamped());
  const auto k = c.clamped();
  CHECK(k.is_clamped());
  CHECK(k.domain_start() == 2.0);
  CHECK(k.domain_end() == 5.0);
  CHECK(max_deviation(c, k) < kExact);
  CHECK(distance(k.control_points().front(), c.evaluate(2.0)) < kExact);
  CHECK(distance(k.control_points().back(), c.evaluate(5.0)) < kExact);

  const auto e = c.elevate_degree(1);
  CHECK(e.degree() == 3);
  CHECK(e.domain_start() == 2.0);
  CHECK(max_deviation(c, e) < kExact);

  const auto same = rational_cubic().clamped();
  CHECK(same.knots() == rational_cubic().knots());
}

TEST_CASE("elevation by zero is the identity", "[elevate]") {
  const auto c = rational_cubic();
  CHECK(c.elevate_degree(0).knots() == c.knots());
}

TEST_CASE("surface degree elevation in u and v", "[elevate][surface]") {
  const auto profile = make_circle({5, 0, 0}, {0, 1, 0}, 1.5);
  const auto torus = revolve(profile, {0, 0, 0}, {0, 0, 1}, 3 * std::numbers::pi / 2);
  const auto eu = torus.elevate_degree_u(1);
  CHECK(eu.degree_u() == 3);
  CHECK(eu.degree_v() == 2);
  CHECK(eu.num_u() == torus.num_u() + 3);  // 3 Bezier spans in u
  CHECK(max_deviation(torus, eu) < kExact);
  const auto ev = torus.elevate_degree_v(2);
  CHECK(ev.degree_v() == 4);
  CHECK(ev.num_v() == torus.num_v() + 2 * 4);  // 4 Bezier spans in v
  CHECK(max_deviation(torus, ev) < kExact);
}

TEST_CASE("ruled surface between curves of different degree", "[elevate][construct]") {
  const auto line = make_line({-1, 0, 0}, {1, 0, 0});
  const auto arc = make_arc({0, 0, 2}, {-1, 0, 0}, {0, 1, 0}, 1.0, 0.0, std::numbers::pi);
  const auto s = ruled(line, arc);
  CHECK(s.degree_u() == 2);
  for (double u = 0.0; u <= 1.0; u += 0.05) {
    CHECK(distance(s.evaluate(u, 0.0), line.evaluate(u)) < kExact);
    CHECK(distance(s.evaluate(u, 1.0), arc.evaluate(u)) < kExact);
  }
}

TEST_CASE("make_compatible with different degrees, knots and clamping", "[elevate][construct]") {
  NurbsCurve a = rational_cubic();
  NurbsCurve b = unclamped_quadratic();
  const NurbsCurve a0 = a, b0 = b;
  make_compatible(a, b);
  CHECK(a.degree() == 3);
  CHECK(b.degree() == 3);
  CHECK(a.knots().size() == b.knots().size());
  for (std::size_t i = 0; i < a.knots().size(); ++i) CHECK(std::abs(a.knots()[i] - b.knots()[i]) < 1e-10);
  // shapes unchanged up to the linear reparametrisation to [0, 1]
  for (double s = 0.0; s <= 1.0; s += 0.01) {
    CHECK(distance(a.evaluate(s), a0.evaluate(s)) < kExact);
    CHECK(distance(b.evaluate(s), b0.evaluate(2.0 + 3.0 * s)) < kExact);
  }
}
