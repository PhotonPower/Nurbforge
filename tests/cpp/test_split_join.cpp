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
#include "nurbforge/nurbs/join.hpp"
#include "nurbforge/nurbs/surface.hpp"

using namespace nf;

namespace {

constexpr double kExact = 1e-12;

// Max distance between a and b on [u0, u1] (same parametrisation).
double max_dev(const NurbsCurve& a, const NurbsCurve& b, double u0, double u1, std::size_t n = 401) {
  double d = 0.0;
  for (std::size_t i = 0; i < n; ++i) {
    const double u = u0 + (u1 - u0) * static_cast<double>(i) / static_cast<double>(n - 1);
    d = std::max(d, distance(a.evaluate(u), b.evaluate(u)));
  }
  return d;
}
double max_dev(const NurbsCurve& part, const NurbsCurve& whole) {
  return max_dev(part, whole, part.domain_start(), part.domain_end());
}

double max_dev(const NurbsSurface& part, const NurbsSurface& whole, std::size_t n = 31) {
  double d = 0.0;
  const auto du = part.domain_u(), dv = part.domain_v();
  for (std::size_t i = 0; i < n; ++i)
    for (std::size_t j = 0; j < n; ++j) {
      const double u = du[0] + (du[1] - du[0]) * static_cast<double>(i) / static_cast<double>(n - 1);
      const double v = dv[0] + (dv[1] - dv[0]) * static_cast<double>(j) / static_cast<double>(n - 1);
      d = std::max(d, distance(part.evaluate(u, v), whole.evaluate(u, v)));
    }
  return d;
}

void check_same_polygon(const NurbsCurve& a, const NurbsCurve& b, double eps = 1e-10) {
  REQUIRE(a.degree() == b.degree());
  REQUIRE(a.knots().size() == b.knots().size());
  REQUIRE(a.num_control_points() == b.num_control_points());
  for (std::size_t i = 0; i < a.knots().size(); ++i) CHECK(std::abs(a.knots()[i] - b.knots()[i]) < eps);
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

NurbsCurve unclamped_quadratic() {
  return {2, {0, 1, 2, 3, 4, 5, 6, 7}, {{0, 0, 0}, {1, 1, 0}, {2, -1, 0}, {3, 1, 1}, {4, 0, 0}}, {}};
}

NurbsCurve scaled_weights(const NurbsCurve& c, double f) {
  std::vector<double> w = c.weights();
  for (double& x : w) x *= f;
  return {c.degree(), c.knots(), c.control_points(), std::move(w)};
}

NurbsSurface torus270() {
  const auto profile = make_circle({5, 0, 0}, {0, 1, 0}, 1.5);
  return revolve(profile, {0, 0, 0}, {0, 0, 1}, 3 * std::numbers::pi / 2);
}

}  // namespace

// ------------------------------------------------------------------ curves

TEST_CASE("split a curve and join it back", "[split][join]") {
  const double u = GENERATE(0.5, 0.3, 0.7, 0.05);  // plain parameter, double knot, simple knot
  const auto c = rational_cubic();
  const auto [l, r] = c.split(u);
  CHECK(l.domain_start() == 0.0);
  CHECK(l.domain_end() == u);
  CHECK(r.domain_start() == u);
  CHECK(r.domain_end() == 1.0);
  CHECK(l.is_clamped());
  CHECK(r.is_clamped());
  CHECK(max_dev(l, c) < kExact);
  CHECK(max_dev(r, c) < kExact);
  CHECK(distance(l.control_points().back(), r.control_points().front()) < kExact);
  check_same_polygon(join(l, r), c);
}

TEST_CASE("split an unclamped curve", "[split]") {
  const auto c = unclamped_quadratic();
  const auto [l, r] = c.split(3.5);
  CHECK(l.domain_start() == 2.0);
  CHECK(r.domain_end() == 5.0);
  CHECK(max_dev(l, c) < kExact);
  CHECK(max_dev(r, c) < kExact);
  check_same_polygon(join(l, r), c.clamped());
}

TEST_CASE("split rejects parameters outside the open domain", "[split]") {
  const auto c = rational_cubic();
  CHECK_THROWS_AS(c.split(0.0), GeometryError);
  CHECK_THROWS_AS(c.split(1.0), GeometryError);
  CHECK_THROWS_AS(c.split(1.2), GeometryError);
}

TEST_CASE("sub_curve", "[split]") {
  const auto c = rational_cubic();
  const auto s = c.sub_curve(0.2, 0.8);
  CHECK(s.domain_start() == 0.2);
  CHECK(s.domain_end() == 0.8);
  CHECK(max_dev(s, c) < kExact);
  check_same_polygon(c.sub_curve(0.0, 1.0), c);
  CHECK(max_dev(c.sub_curve(0.0, 0.3), c) < kExact);
  CHECK(max_dev(c.sub_curve(0.3, 1.0), c) < kExact);
  CHECK_THROWS_AS(c.sub_curve(0.6, 0.4), GeometryError);
  CHECK_THROWS_AS(c.sub_curve(-0.1, 0.4), GeometryError);
}

TEST_CASE("join curves of different degree", "[join]") {
  const auto line = make_line({-1, 0, 0}, {1, 0, 0});
  const auto arc = make_arc({0, 0, 0}, {1, 0, 0}, {0, 1, 0}, 1.0, 0.0, std::numbers::pi);
  const auto j = join(line, arc);
  CHECK(j.degree() == 2);
  CHECK(j.domain_start() == 0.0);
  CHECK(j.domain_end() == 2.0);
  CHECK(knot_multiplicity(j.knots(), 1.0) == 2);  // corner stays a corner
  CHECK(max_dev(j, line, 0.0, 1.0) < kExact);
  for (double u = 0.0; u <= 1.0; u += 0.01) CHECK(distance(j.evaluate(1.0 + u), arc.evaluate(u)) < kExact);
}

TEST_CASE("join rescales weights of rational curves", "[join]") {
  const auto circle = make_circle({1, 1, 0}, {0, 0, 1}, 2.0);
  const auto [l, r] = circle.split(0.4);
  const auto j = join(scaled_weights(l, 0.25), scaled_weights(r, 7.0));
  CHECK(max_dev(j, circle, 0.0, 1.0) < kExact);
  CHECK(j.num_control_points() == circle.num_control_points());  // 0.4 removed again
}

TEST_CASE("join closes small gaps and rejects large ones", "[join]") {
  const auto a = make_line({0, 0, 0}, {1, 0, 0});
  const auto b = make_line({1, 1e-6, 0}, {2, 0, 0});
  CHECK_THROWS_AS(join(a, b), GeometryError);
  const auto j = join(a, b, 1e-5);
  CHECK(distance(j.evaluate(1.0), {1, 0.5e-6, 0}) < kExact);
  CHECK(distance(j.evaluate(0.0), {0, 0, 0}) < kExact);
  CHECK(distance(j.evaluate(2.0), {2, 0, 0}) < kExact);
}

TEST_CASE("joining collinear lines removes the joint knot", "[join]") {
  const auto j = join(make_line({0, 0, 0}, {1, 1, 1}), make_line({1, 1, 1}, {2, 2, 2}));
  CHECK(j.num_control_points() == 2);
  CHECK(distance(j.evaluate(1.0), {1, 1, 1}) < kExact);
}

// ------------------------------------------------------------------ surfaces

TEST_CASE("split a surface in u and v and join it back", "[split][join][surface]") {
  const auto t = torus270();
  const auto [a, b] = t.split_u(0.4);
  CHECK(a.domain_u()[1] == 0.4);
  CHECK(b.domain_u()[0] == 0.4);
  CHECK(max_dev(a, t) < kExact);
  CHECK(max_dev(b, t) < kExact);
  const auto ju = join_u(a, b);
  CHECK(ju.num_u() == t.num_u());
  CHECK(max_dev(t, ju) < kExact);

  const auto [c, d] = t.split_v(0.6);
  CHECK(max_dev(c, t) < kExact);
  CHECK(max_dev(d, t) < kExact);
  const auto jv = join_v(c, d);
  CHECK(jv.num_v() == t.num_v());
  CHECK(max_dev(t, jv) < kExact);
}

TEST_CASE("sub_surface and clamped surface", "[split][surface]") {
  const auto t = torus270();
  const auto s = t.sub_surface(0.1, 0.55, 0.2, 0.9);
  CHECK(s.domain_u()[0] == 0.1);
  CHECK(s.domain_v()[1] == 0.9);
  CHECK(max_dev(s, t) < kExact);
  CHECK_THROWS_AS(t.sub_surface(0.5, 0.4, 0.0, 1.0), GeometryError);

  const auto e = extrude(unclamped_quadratic(), {0, 0, 3});
  const auto ec = e.clamped();
  CHECK(ec.knots_u().front() == ec.domain_u()[0]);
  CHECK(max_dev(e, ec) < kExact);
}

TEST_CASE("join_u makes the v directions compatible", "[join][surface]") {
  const auto a = make_rectangle({0, 0, 0}, {1, 0, 0}, {0, 1, 0});
  // b: adjacent patch with a different v structure (degree 2, extra knot, other v domain)
  const auto b0 = make_rectangle({1, 0, 0}, {1, 0, 0.5}, {0, 1, 0}).elevate_degree_v(1).refine_knots_v({0.3});
  const KnotVector kv{5, 5, 5, 5.9, 8, 8, 8};  // same knot 0.3, mapped to the v domain [5, 8]
  const NurbsSurface b{b0.degree_u(), b0.degree_v(), b0.knots_u(), kv, b0.num_u(), b0.num_v(),
                       b0.control_points(), b0.weights()};
  const auto j = join_u(a, b);
  CHECK(j.domain_u()[1] == 2.0);
  CHECK(j.degree_v() == 2);
  for (double u = 0.0; u <= 1.0; u += 0.1)
    for (double v = 0.0; v <= 1.0; v += 0.1) {
      CHECK(distance(j.evaluate(u, v), a.evaluate(u, v)) < kExact);
      CHECK(distance(j.evaluate(1.0 + u, v), b.evaluate(u, 5.0 + 3.0 * v)) < kExact);
    }
  CHECK_THROWS_AS(join_u(a, a), GeometryError);  // boundaries far apart
}
