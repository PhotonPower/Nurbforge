#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <numbers>
#include <vector>

#include "nurbforge/analysis/closest_point.hpp"
#include "nurbforge/core/error.hpp"
#include "nurbforge/geometry/primitives.hpp"
#include "nurbforge/geometry/transform.hpp"
#include "nurbforge/modeling/construct.hpp"
#include "nurbforge/modeling/loft.hpp"

using namespace nf;
constexpr double pi = std::numbers::pi;

namespace {

std::vector<NurbsCurve> circles(const std::vector<double>& z, const std::vector<double>& r) {
  std::vector<NurbsCurve> s;
  for (std::size_t k = 0; k < z.size(); ++k) s.push_back(make_circle({0, 0, z[k]}, {0, 0, 1}, r[k]));
  return s;
}

// Every section lies on the surface at its parameter: S(u, v_k) = C_k(u) after making the
// sections compatible (domain [0, 1]).
void check_sections(const NurbsSurface& s, const std::vector<NurbsCurve>& sections,
                    const std::vector<double>& v) {
  std::vector<NurbsCurve> c = sections;
  make_compatible(c);
  for (std::size_t k = 0; k < c.size(); ++k)
    for (double u = 0.0; u <= 1.0; u += 0.05) CHECK(distance(s.evaluate(u, v[k]), c[k].evaluate(u)) < 1e-11);
}

}  // namespace

TEST_CASE("loft through equal circles is an exact cylinder", "[loft]") {
  const auto sec = circles({0, 1, 2.5, 4}, {2, 2, 2, 2});
  const auto s = loft(sec);
  CHECK(s.degree_v() == 3);
  check_sections(s, sec, loft_parameters(sec));
  for (double u = 0.0; u <= 1.0; u += 0.05)
    for (double v = 0.0; v <= 1.0; v += 0.05) {
      const Vec3 p = s.evaluate(u, v);
      CHECK(std::abs(std::hypot(p.x, p.y) - 2.0) < 1e-12);
    }
}

TEST_CASE("loft through linearly growing circles is an exact cone", "[loft]") {
  const auto sec = circles({0, 1, 2, 3, 4}, {1, 1.5, 2, 2.5, 3});
  const auto s = loft(sec, 3, std::vector<double>{0, 0.25, 0.5, 0.75, 1});
  for (double u = 0.0; u <= 1.0; u += 0.05)
    for (double v = 0.0; v <= 1.0; v += 0.05) {
      const Vec3 p = s.evaluate(u, v);
      CHECK(std::abs(std::hypot(p.x, p.y) - (1.0 + 0.5 * p.z)) < 1e-12);
    }
}

TEST_CASE("loft through sections of different degree and knots", "[loft]") {
  const auto a = NurbsCurve::from_points(2, {{0, 0, 0}, {1, 1, 0}, {2, 0, 0}, {3, 1, 0}});
  const auto b = make_line({0, 0, 1}, {3, 0, 1});
  const auto c = NurbsCurve::from_points(3, {{0, 0, 2}, {1, -1, 2}, {1.5, 2, 2}, {2, -1, 2}, {3, 0, 2}});
  const std::vector<NurbsCurve> sec{a, b, c};
  const auto s = loft(sec, 2);
  CHECK(s.degree_u() == 3);
  CHECK(s.degree_v() == 2);
  check_sections(s, sec, loft_parameters(sec));
}

TEST_CASE("loft through latitude circles approximates a sphere", "[loft]") {
  std::vector<NurbsCurve> sec;
  for (int k = 1; k < 10; ++k) {
    const double th = pi * k / 10.0;
    sec.push_back(make_circle({0, 0, -std::cos(th)}, {0, 0, 1}, std::sin(th)));
  }
  const auto s = loft(sec, 3, Parametrization::ChordLength);
  check_sections(s, sec, loft_parameters(sec, Parametrization::ChordLength));
  for (double u = 0.0; u <= 1.0; u += 0.1)
    for (double v = 0.0; v <= 1.0; v += 0.05) CHECK(std::abs(norm(s.evaluate(u, v)) - 1.0) < 2e-3);
}

TEST_CASE("linear loft through two sections equals the ruled surface", "[loft]") {
  const auto line = make_line({-1, 0, 0}, {1, 0, 0});
  const auto arc = make_arc({0, 0, 2}, {-1, 0, 0}, {0, 1, 0}, 1.0, 0.0, pi);
  const auto l = loft({line, arc}, 1);
  const auto r = ruled(line, arc);
  for (double u = 0.0; u <= 1.0; u += 0.1)
    for (double v = 0.0; v <= 1.0; v += 0.1) CHECK(distance(l.evaluate(u, v), r.evaluate(u, v)) < 1e-12);
}

TEST_CASE("loft of a twisted stack of rotated profiles", "[loft]") {
  const auto profile = make_polyline({{1, 0, 0}, {2, 0.5, 0}, {1.5, 2, 0}, {0.5, 1.5, 0}});
  std::vector<NurbsCurve> sec;
  for (int k = 0; k < 6; ++k)
    sec.push_back(profile.transformed(Transform::translation({0, 0, double(k)}) *
                                      Transform::rotation({0, 0, 0}, {0, 0, 1}, 0.2 * k)));
  const auto s = loft(sec);
  check_sections(s, sec, loft_parameters(sec));
  // a point of the surface between two sections lies between their heights
  const Vec3 p = s.evaluate(0.5, 0.5);
  CHECK(p.z > 2.0);
  CHECK(p.z < 3.0);
}

TEST_CASE("loft input checks", "[loft]") {
  const auto c = make_circle({0, 0, 0}, {0, 0, 1}, 1.0);
  CHECK_THROWS_AS(loft({c}), GeometryError);
  CHECK_THROWS_AS(loft({c, c}), GeometryError);  // coincident sections
  CHECK_THROWS_AS(loft({c, c.transformed(Transform::translation({0, 0, 1}))}, 3, std::vector<double>{0, 0}),
                  GeometryError);
}
