#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <numbers>

#include "nurbforge/analysis/properties.hpp"
#include "nurbforge/analysis/section.hpp"
#include "nurbforge/core/error.hpp"
#include "nurbforge/geometry/primitives.hpp"
#include "nurbforge/modeling/construct.hpp"

using namespace nf;
constexpr double pi = std::numbers::pi;

TEST_CASE("tilted planar section of a sphere", "[section]") {
  const double R = 2.0, d = 0.9;
  const auto sphere = revolve(make_arc({0, 0, 0}, {0, 0, -1}, {1, 0, 0}, R, 0.0, pi), {0, 0, 0}, {0, 0, 1}, 2 * pi);
  const Vec3 n = normalized(Vec3{1, 2, 2});
  const Plane plane{d * n, n};  // plane.point is the foot of the centre: circle about (0, 0)
  const auto r = planar_section(sphere, plane);
  const double rho = std::sqrt(R * R - d * d);
  double L = 0.0;
  for (const auto& c : r.curves) {
    L += arc_length(c.curve);
    for (const auto& p : c.points) {
      CHECK(std::abs(plane.distance(p.point)) < 1e-9);
      CHECK(std::abs(std::hypot(p.u2, p.v2) - rho) < 1e-9);  // in-plane coordinates
    }
    for (double f = 0.0; f <= 1.0; f += 0.1) {
      const double t = c.curve.domain_start() + f * (c.curve.domain_end() - c.curve.domain_start());
      const Vec3 st = c.pcurve2.evaluate(t);
      CHECK(std::abs(std::hypot(st.x, st.y) - rho) < 1e-6);
    }
  }
  CHECK(std::abs(L - 2 * pi * rho) < 1e-5);
}

TEST_CASE("parallel sections of a torus", "[section]") {
  const double R = 5.0, rr = 1.5;
  const auto torus = revolve(make_circle({R, 0, 0}, {0, 1, 0}, rr), {0, 0, 0}, {0, 0, 1}, 2 * pi);
  const auto slices = planar_sections(torus, {{0, 0, -1.2}, {0, 0, 1}}, 0.6, 5);  // z = -1.2 .. 1.2
  REQUIRE(slices.size() == 5);
  for (std::size_t k = 0; k < slices.size(); ++k) {
    const double z = -1.2 + 0.6 * static_cast<double>(k);
    const double w = std::sqrt(rr * rr - z * z);
    double L = 0.0;
    for (const auto& c : slices[k].curves) {
      L += arc_length(c.curve);
      for (const auto& p : c.points) CHECK(std::abs(p.point.z - z) < 1e-9);
    }
    CHECK(std::abs(L - 2 * pi * (R + w) - 2 * pi * (R - w)) < 1e-4);
  }
  CHECK_THROWS_AS(planar_sections(torus, {{0, 0, 0}, {0, 0, 1}}, 0.0, 3), GeometryError);
}

TEST_CASE("curve / plane intersection", "[section]") {
  const auto circle = make_circle({1, 0, 0}, {1, 0, 1}, 2.0);
  const Plane plane{{0, 0, 0.5}, {0, 0, 1}};
  const auto hits = intersect(circle, plane);
  REQUIRE(hits.size() == 2);
  for (const auto& h : hits) {
    CHECK(std::abs(h.point.z - 0.5) < 1e-12);
    CHECK(std::abs(distance(h.point, {1, 0, 0}) - 2.0) < 1e-12);
  }
  CHECK(intersect(make_line({0, 0, 1}, {1, 1, 1}), plane).empty());
}

TEST_CASE("curve normal to the plane", "[section]") {
  const auto hits = intersect(make_line({0, 0, -3}, {0, 0, 3}), Plane{{0, 0, 0.5}, {0, 0, 1}});
  REQUIRE(hits.size() == 1);
  CHECK(std::abs(hits[0].point.z - 0.5) < 1e-12);
  CHECK(std::abs(hits[0].u) < 1e-12);  // in-plane coordinates of the foot (0, 0)
  CHECK(std::abs(hits[0].v) < 1e-12);
}
