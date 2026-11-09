#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <numbers>
#include <vector>

#include "nurbforge/analysis/properties.hpp"
#include "nurbforge/core/error.hpp"
#include "nurbforge/geometry/primitives.hpp"
#include "nurbforge/modeling/fillet.hpp"
#include "nurbforge/nurbs/knot_vector.hpp"

using namespace nf;
using Catch::Approx;
constexpr double pi = std::numbers::pi;

namespace {

double signed_area_xy(const NurbsCurve& c) {  // Green, Gauss-Legendre on fine sub-intervals
  static const double gx[5] = {-0.9061798459386640, -0.5384693101056831, 0.0, 0.5384693101056831,
                               0.9061798459386640};
  static const double gw[5] = {0.2369268850561891, 0.4786286704993665, 0.5688888888888889,
                               0.4786286704993665, 0.2369268850561891};
  const auto k = unique_knots(c.knots()).first;
  double a = 0;
  for (std::size_t i = 0; i + 1 < k.size(); ++i)
    for (int q = 0; q < 8; ++q) {
      const double a0 = k[i] + (k[i + 1] - k[i]) * q / 8,
                   a1 = k[i] + (k[i + 1] - k[i]) * (q + 1) / 8;
      for (int g = 0; g < 5; ++g) {
        const auto d = c.derivatives(0.5 * (a0 + a1) + 0.5 * (a1 - a0) * gx[g], 1);
        a += 0.25 * (a1 - a0) * gw[g] * (d[0].x * d[1].y - d[0].y * d[1].x);
      }
    }
  return a;
}

/// Largest angle between the one-sided tangents at the interior knots.
double max_kink(const NurbsCurve& c) {
  const auto k = unique_knots(c.knots()).first;
  double m = 0;
  for (std::size_t i = 1; i + 1 < k.size(); ++i) {
    const Vec3 l = normalized(c.sub_curve(k[i - 1], k[i]).derivatives(k[i], 1)[1]);
    const Vec3 r = normalized(c.sub_curve(k[i], k[i + 1]).derivatives(k[i], 1)[1]);
    m = std::max(m, std::atan2(norm(cross(l, r)), dot(l, r)));
  }
  return m;
}

}  // namespace

TEST_CASE("circle and arc through 3 points", "[constructors]") {
  const Vec3 p0{1, 0, 2}, p1{0, 3, 1}, p2{-2, 1, 0};
  const auto c = make_circle_3pt(p0, p1, p2);
  CHECK(distance(c.evaluate(c.domain_start()), p0) < 1e-12);
  // all points on one circle: equal distance from the centre of three samples' circumcircle
  const auto s = c.sample(50);
  double rmin = 1e300, rmax = 0;
  // centre: the mean of opposite points
  const Vec3 centre =
      0.5 * (c.evaluate(0.0) + c.evaluate(0.5 * (c.domain_start() + c.domain_end())));
  for (const auto& p : s) {
    rmin = std::min(rmin, distance(p, centre));
    rmax = std::max(rmax, distance(p, centre));
  }
  CHECK(rmax - rmin < 1e-12);
  CHECK(std::abs(distance(p1, centre) - rmax) < 1e-12);
  CHECK(std::abs(distance(p2, centre) - rmax) < 1e-12);

  const auto arc = make_arc_3pt(p0, p1, p2);
  CHECK(distance(arc.evaluate(arc.domain_start()), p0) == 0.0);
  CHECK(distance(arc.evaluate(arc.domain_end()), p2) == 0.0);
  // p1 lies on the arc
  double best = 1e300;
  for (const auto& p : arc.sample(2001)) best = std::min(best, distance(p, p1));
  CHECK(best < 1e-2);
  // a major arc (> 180 degrees)
  const auto major = make_arc_3pt({1, 0, 0}, {-1, 0, 0}, {0, -1, 0});
  CHECK(arc_length(major) == Approx(1.5 * pi).epsilon(1e-10));
  CHECK_THROWS_AS(make_arc_3pt({0, 0, 0}, {1, 1, 1}, {2, 2, 2}), GeometryError);
}

TEST_CASE("polygon and rectangle outlines", "[constructors]") {
  const auto hex = make_polygon({1, 1, 0}, {0, 0, 1}, 2.0, 6, {1, 0, 0});
  CHECK(arc_length(hex) == Approx(12.0).epsilon(1e-12));
  CHECK(distance(hex.evaluate(hex.domain_start()), Vec3{3, 1, 0}) < 1e-14);
  CHECK(signed_area_xy(hex) == Approx(1.5 * std::sqrt(3.0) * 4).epsilon(1e-12));
  const auto r = make_rectangle_curve({0, 0, 0}, {3, 0, 0}, {0, 2, 0});
  CHECK(signed_area_xy(r) == Approx(6.0).epsilon(1e-12));
  CHECK_THROWS_AS(make_polygon({0, 0, 0}, {0, 0, 1}, 1.0, 2), GeometryError);
}

TEST_CASE("helix and spiral", "[constructors]") {
  SECTION("cylindrical, right-handed") {
    const double r = 2, pitch = 0.5, turns = 3.25;
    const auto h = make_helix({0, 0, 0}, {0, 0, 1}, {r, 0, 1}, pitch, turns, -1, false, 1e-8);
    CHECK(h.degree() == 5);
    CHECK(h.domain_end() == turns);
    for (int i = 0; i <= 400; ++i) {
      const double t = turns * i / 400.0, a = 2 * pi * t;
      const Vec3 exact{r * std::cos(a), r * std::sin(a), 1 + pitch * t};
      CHECK(distance(h.evaluate(t), exact) < 1e-8);
    }
    CHECK(arc_length(h) == Approx(turns * std::hypot(2 * pi * r, pitch)).epsilon(1e-8));
    // C2: one-sided second derivatives agree at the interior knots
    const auto k = unique_knots(h.knots());
    for (std::size_t i = 1; i + 1 < k.first.size(); ++i) CHECK(k.second[i] == 3);
  }
  SECTION("conical, left-handed, oblique axis") {
    const Vec3 a = normalized(Vec3{1, 2, 2}), o{1, -1, 0.5};
    const Vec3 x = any_perpendicular(a);
    const auto h = make_helix(o, a, o + 3.0 * x, 1.0, 2.0, 1.0, true, 1e-8);
    const Vec3 y = cross(x, a);  // clockwise about a
    for (int i = 0; i <= 200; ++i) {
      const double t = 2.0 * i / 200.0, ang = 2 * pi * t, rr = 3.0 - t;
      const Vec3 exact = o + t * a + rr * (std::cos(ang) * x + std::sin(ang) * y);
      CHECK(distance(h.evaluate(t), exact) < 1e-8);
    }
  }
  SECTION("planar Archimedean spiral from the axis") {
    const auto s = make_helix({0, 0, 0}, {0, 0, 1}, {0, 0, 0}, 0.0, 4.0, 4.0);
    for (const auto& p : s.sample(100)) CHECK(std::abs(p.z) < 1e-12);
    CHECK(norm(s.evaluate(4.0)) == Approx(4.0).epsilon(1e-9));
    CHECK(norm(s.evaluate(2.5)) == Approx(2.5).epsilon(1e-6));
  }
  CHECK_THROWS_AS(make_helix({0, 0, 0}, {0, 0, 1}, {1, 0, 0}, 1, 0), GeometryError);
}

TEST_CASE("fillet_corners", "[constructors]") {
  const double r = 0.3;
  SECTION("closed square: rounded rectangle") {
    const auto sq = make_rectangle_curve({0, 0, 0}, {2, 0, 0}, {0, 2, 0});
    const auto c = fillet_corners(sq, r);
    CHECK(distance(c.evaluate(c.domain_start()), c.evaluate(c.domain_end())) < 1e-9);
    CHECK(arc_length(c) == Approx(8 - 8 * r + 2 * pi * r).epsilon(1e-10));
    CHECK(signed_area_xy(c) == Approx(4 - (4 - pi) * r * r).epsilon(1e-10));
    CHECK(max_kink(c) < 1e-8);
  }
  SECTION("open L-polyline with a reflex and a convex corner") {
    const auto pl = make_polyline({{0, 0, 0}, {2, 0, 0}, {2, 1, 0}, {3, 1, 0}});
    const auto c = fillet_corners(pl, r);
    CHECK(distance(c.evaluate(c.domain_start()), Vec3{0, 0, 0}) < 1e-12);
    CHECK(distance(c.evaluate(c.domain_end()), Vec3{3, 1, 0}) < 1e-12);
    CHECK(arc_length(c) == Approx(4 - 4 * r + pi * r).epsilon(1e-10));
    CHECK(max_kink(c) < 1e-8);
  }
  SECTION("hexagon in a tilted plane") {
    const auto hex = make_polygon({0, 0, 0}, {0, 1, 1}, 2.0, 6);
    const auto c = fillet_corners(hex, r);
    // each 60-degree corner cuts 2 r tan(30) and adds the arc r pi / 3
    CHECK(arc_length(c) == Approx(12 - 12 * r * std::tan(pi / 6) + 2 * pi * r).epsilon(1e-10));
    CHECK(max_kink(c) < 1e-8);
  }
  SECTION("errors and no-ops") {
    CHECK_THROWS_AS(fillet_corners(make_rectangle_curve({0, 0, 0}, {1, 0, 0}, {0, 1, 0}), 0.6),
                    GeometryError);
    CHECK_THROWS_AS(
        fillet_corners(make_polyline({{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {1, 1, 1}}), 0.1),
        GeometryError);  // not planar
    const auto line = make_line({0, 0, 0}, {1, 1, 0});
    CHECK(fillet_corners(line, 0.2).control_points().size() == 2);
  }
}
