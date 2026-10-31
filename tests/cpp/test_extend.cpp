#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <algorithm>
#include <cmath>
#include <numbers>

#include "nurbforge/analysis/properties.hpp"
#include "nurbforge/core/error.hpp"
#include "nurbforge/geometry/primitives.hpp"
#include "nurbforge/modeling/construct.hpp"
#include "nurbforge/modeling/extend.hpp"

using namespace nf;
constexpr double pi = std::numbers::pi;

namespace {

NurbsCurve wavy_cubic() {
  return NurbsCurve::from_points(
      3, {{0, 0, 0}, {1, 3, 0}, {2, -2, 0}, {3, 3, 1}, {4, -1, 0}, {5, 2, -1}, {6, 0, 0}});
}

NurbsCurve rational_cubic() {
  return {3,
          {0, 0, 0, 0, 0.3, 0.3, 0.7, 1, 1, 1, 1},
          {{0, 0, 0}, {1, 2, 0}, {2, 3, 1}, {3, 1, 1}, {4, 0, 2}, {5, 2, 2}, {6, 1, 0}},
          {1.0, 0.5, 2.0, 1.0, 0.7, 1.3, 1.0}};
}

// The original part of `e` coincides with `c` at the same parameters.
double original_deviation(const NurbsCurve& c, const NurbsCurve& e) {
  double d = 0.0;
  for (int i = 0; i <= 200; ++i) {
    const double u = c.domain_start() + (c.domain_end() - c.domain_start()) * i / 200.0;
    d = std::max(d, distance(c.evaluate(u), e.evaluate(u)));
  }
  return d;
}

}  // namespace

TEST_CASE("natural curve extension: smooth continuation of the end piece", "[extend]") {
  const auto c = GENERATE(wavy_cubic(), rational_cubic());
  const double L = 0.8;
  const auto e = extend(c, CurveEnd::End, L);
  CHECK(e.domain_start() == c.domain_start());
  CHECK(e.domain_end() > c.domain_end());
  CHECK(original_deviation(c, e) < 1e-12);
  CHECK(std::abs(arc_length(e, c.domain_end(), e.domain_end()) - L) < 1e-9);
  // C-infinity at the old end: derivatives of the old end piece agree up to order 3
  const double b = c.domain_end();
  const auto dc = c.derivatives(b, 3), de = e.derivatives(b, 3);
  for (std::size_t k = 1; k <= 3; ++k) CHECK(distance(dc[k], de[k]) < 1e-8 * (1 + norm(dc[k])));
  CHECK(e.num_control_points() == c.num_control_points());  // only the end piece changed
}

TEST_CASE("natural extension at the start keeps parameters", "[extend]") {
  const auto c = wavy_cubic();
  const auto e = extend(c, CurveEnd::Start, 1.3);
  CHECK(e.domain_end() == c.domain_end());
  CHECK(e.domain_start() < 0.0);
  CHECK(original_deviation(c, e) < 1e-12);
  CHECK(std::abs(arc_length(e, e.domain_start(), 0.0) - 1.3) < 1e-9);
}

TEST_CASE("natural extension of a circular arc stays on the circle", "[extend]") {
  const double R = 2.0;
  const auto arc = make_arc({1, 1, 0}, {1, 0, 0}, {0, 1, 0}, R, 0.0, pi / 2);
  const auto e = extend(extend(arc, CurveEnd::End, R * pi / 6), CurveEnd::Start, R * pi / 8);
  for (int i = 0; i <= 100; ++i) {
    const double u = e.domain_start() + (e.domain_end() - e.domain_start()) * i / 100.0;
    CHECK(std::abs(distance(e.evaluate(u), {1, 1, 0}) - R) < 1e-12);
  }
  CHECK(std::abs(arc_length(e) - R * (pi / 2 + pi / 6 + pi / 8)) < 1e-9);
  // continuing a quarter arc of 90 degrees by another 120 degrees needs non-positive weights
  CHECK_THROWS_AS(extend(arc, CurveEnd::End, R * 2 * pi / 3), GeometryError);
}

TEST_CASE("linear curve extension", "[extend]") {
  const auto c = wavy_cubic();
  const double L = 2.0;
  const auto e = extend(c, CurveEnd::End, L, Extension::Linear);
  const double b = c.domain_end();
  const Vec3 P = c.evaluate(b), T = c.tangent(b);
  CHECK(distance(e.evaluate(e.domain_end()), P + L * T) < 1e-12);
  CHECK(original_deviation(c, e) < 1e-12);
  for (double s = 0.1; s < 1.0; s += 0.2) {  // straight
    const Vec3 q = e.evaluate(b + s * (e.domain_end() - b));
    CHECK(norm(cross(q - P, T)) < 1e-12);
  }
  // C1 for non-rational curves
  CHECK(distance(c.derivatives(b, 1)[1], e.derivatives(b, 1)[1]) < 1e-10);
  const auto s = extend(c, CurveEnd::Start, L, Extension::Linear);
  CHECK(distance(s.evaluate(s.domain_start()), c.evaluate(0.0) - L * c.tangent(0.0)) < 1e-12);
  CHECK(original_deviation(c, s) < 1e-12);
}

TEST_CASE("extending a line", "[extend]") {
  const auto l = make_line({0, 0, 0}, {3, 4, 0});  // length 5, speed 5
  for (auto type : {Extension::Natural, Extension::Linear}) {
    const auto e = extend(l, CurveEnd::End, 5.0, type);
    CHECK(std::abs(e.domain_end() - 2.0) < 1e-12);
    CHECK(distance(e.evaluate(e.domain_end()), {6, 8, 0}) < 1e-12);
  }
  CHECK_THROWS_AS(extend(l, CurveEnd::End, 0.0), GeometryError);
}

TEST_CASE("natural surface extension of a torus stays on the torus", "[extend][surface]") {
  const double R = 5.0, r = 1.5;
  const auto torus = revolve(make_circle({R, 0, 0}, {0, 1, 0}, r), {0, 0, 0}, {0, 0, 1}, pi);
  for (auto side : {SurfaceSide::UMax, SurfaceSide::UMin, SurfaceSide::VMax, SurfaceSide::VMin}) {
    // u: around the axis (half turn), v: around the tube (full circle, extend a bit)
    const auto e = extend(torus, side, 0.4);
    const auto du = e.domain_u(), dv = e.domain_v();
    for (int i = 0; i <= 20; ++i)
      for (int j = 0; j <= 20; ++j) {
        const Vec3 p = e.evaluate(du[0] + (du[1] - du[0]) * i / 20.0, dv[0] + (dv[1] - dv[0]) * j / 20.0);
        const double rho = std::hypot(p.x, p.y);
        CHECK(std::abs(std::hypot(rho - R, p.z) - r) < 1e-10);
      }
    // original parameters unchanged
    CHECK(distance(e.evaluate(0.3, 0.6), torus.evaluate(0.3, 0.6)) < 1e-12);
  }
}

TEST_CASE("surface extension along the extrusion direction", "[extend][surface]") {
  const auto cyl = extrude(make_circle({0, 0, 0}, {0, 0, 1}, 2.0), {0, 0, 4});
  const auto e = extend(cyl, SurfaceSide::VMax, 1.0);
  CHECK(std::abs(e.evaluate(0.3, e.domain_v()[1]).z - 5.0) < 1e-12);
  const auto f = extend(cyl, SurfaceSide::VMin, 1.5, Extension::Linear);
  CHECK(std::abs(f.evaluate(0.7, f.domain_v()[0]).z + 1.5) < 1e-12);
  CHECK_THROWS_AS(extend(cyl, SurfaceSide::UMax, 1.0, Extension::Linear), GeometryError);  // rational
}

TEST_CASE("linear surface extension is C1", "[extend][surface]") {
  std::vector<Vec3> P;
  for (int i = 0; i < 5; ++i)
    for (int j = 0; j < 4; ++j) P.push_back({double(i), double(j), std::sin(1.1 * i) * std::cos(0.9 * j)});
  const NurbsSurface s{3, 3, clamped_uniform_knots(3, 5), clamped_uniform_knots(3, 4), 5, 4, P};
  const auto e = extend(s, SurfaceSide::UMax, 1.0, Extension::Linear);
  for (double v = 0.0; v <= 1.0; v += 0.1) {
    CHECK(distance(e.evaluate(1.0, v), s.evaluate(1.0, v)) < 1e-12);
    CHECK(distance(e.derivatives(1.0, v, 1)[1][0], s.derivatives(1.0, v, 1)[1][0]) < 1e-10);
  }
  CHECK(distance(e.evaluate(0.4, 0.7), s.evaluate(0.4, 0.7)) < 1e-12);
  // mid iso-curve extended by the requested length
  const auto iso = e.iso_curve_v(0.5);
  CHECK(std::abs(arc_length(iso, 1.0, e.domain_u()[1]) - 1.0) < 1e-9);
}
