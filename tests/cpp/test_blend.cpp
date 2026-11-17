#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <numbers>

#include "nurbforge/analysis/properties.hpp"
#include "nurbforge/core/error.hpp"
#include "nurbforge/geometry/primitives.hpp"
#include "nurbforge/modeling/blend.hpp"

using namespace nf;
constexpr double pi = std::numbers::pi;

namespace {

Vec3 unit_tangent(const NurbsCurve& c, double u) { return normalized(c.derivatives(u, 1)[1]); }
Vec3 curvature_vector(const NurbsCurve& c, double u) {
  const auto d = c.derivatives(u, 2);
  const double s = norm(d[1]);
  return (1.0 / (s * s * s * s)) * cross(cross(d[1], d[2]), d[1]);
}

}  // namespace

TEST_CASE("blend curves: G1 and G2 transitions", "[blend]") {
  SECTION("two lines, G1: tangents continue") {
    const auto a = make_line({0, 0, 0}, {1, 0, 0}), b = make_line({2, 1, 0}, {2, 2, 0});
    const auto bl = blend_curves(a, CurveEnd::End, b, CurveEnd::Start, 1);
    CHECK(bl.degree() == 3);
    CHECK(distance(bl.evaluate(0.0), Vec3{1, 0, 0}) <= 1e-15);
    CHECK(distance(bl.evaluate(1.0), Vec3{2, 1, 0}) <= 1e-15);
    CHECK(distance(unit_tangent(bl, 0.0), Vec3{1, 0, 0}) <= 1e-14);
    CHECK(distance(unit_tangent(bl, 1.0), Vec3{0, 1, 0}) <= 1e-14);
  }
  SECTION("two arcs, G2: tangents and curvature vectors continue") {
    const double R = 2.0;
    const auto a = make_arc({0, 0, 0}, {1, 0, 0}, {0, 1, 0}, R, 0.0, pi / 4);
    const auto b = make_arc({5, 0, 0}, {1, 0, 0}, {0, 1, 0}, 1.0, pi / 2, pi);  // left half
    for (const auto e1 : {CurveEnd::Start, CurveEnd::End}) {
      const auto bl = blend_curves(a, CurveEnd::End, b, e1, 2);
      CHECK(bl.degree() == 5);
      const double ub = e1 == CurveEnd::Start ? b.domain_start() : b.domain_end();
      const double sign = e1 == CurveEnd::Start ? 1.0 : -1.0;
      CHECK(distance(unit_tangent(bl, 0.0), unit_tangent(a, a.domain_end())) <= 1e-12);
      CHECK(distance(unit_tangent(bl, 1.0), sign * unit_tangent(b, ub)) <= 1e-12);
      CHECK(distance(curvature_vector(bl, 0.0), curvature_vector(a, a.domain_end())) <= 1e-9);
      CHECK(distance(curvature_vector(bl, 1.0), curvature_vector(b, ub)) <= 1e-9);
      CHECK(std::abs(curvature(bl, 0.0) - 1.0 / R) <= 1e-9);
    }
  }
  SECTION("leaving a curve at its start; argument checks") {
    const auto a = make_line({0, 0, 0}, {1, 0, 0}), b = make_line({-1, 1, 0}, {-2, 1, 0});
    const auto bl = blend_curves(a, CurveEnd::Start, b, CurveEnd::Start, 2);
    CHECK(distance(unit_tangent(bl, 0.0), Vec3{-1, 0, 0}) <= 1e-14);
    CHECK(norm(curvature_vector(bl, 0.0)) <= 1e-12);  // straight ends: zero curvature (G2)
    CHECK_THROWS_AS(blend_curves(a, CurveEnd::End, a, CurveEnd::End, 2), GeometryError);
    CHECK_THROWS_AS(blend_curves(a, CurveEnd::End, b, CurveEnd::Start, 3), GeometryError);
  }
}
