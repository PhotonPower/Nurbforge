#include <algorithm>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <numbers>
#include <utility>
#include <vector>

#include "nurbforge/analysis/closest_point.hpp"
#include "nurbforge/analysis/hatch.hpp"
#include "nurbforge/analysis/properties.hpp"
#include "nurbforge/core/error.hpp"
#include "nurbforge/geometry/primitives.hpp"
#include "nurbforge/geometry/transform.hpp"
#include "nurbforge/modeling/construct.hpp"
#include "nurbforge/topology/shell.hpp"
#include "nurbforge/topology/trimmed_surface.hpp"

using namespace nf;
using Catch::Approx;
constexpr double pi = std::numbers::pi;

namespace {

using Loop = std::vector<NurbsCurve>;

Loop circle(double cx, double cy, double r, bool ccw = true) {
  const auto c = make_circle({cx, cy, 0}, {0, 0, 1}, r);
  return {ccw ? c : c.reversed()};
}

Loop polygon(const std::vector<Vec3>& pts) {  // closed, in the given order
  Loop l;
  for (std::size_t i = 0; i < pts.size(); ++i)
    l.push_back(make_line(pts[i], pts[(i + 1) % pts.size()]));
  return l;
}

Loop square(double cx, double cy, double h, double angle = 0.0, bool ccw = true) {
  std::vector<Vec3> p;
  for (int i = 0; i < 4; ++i) {
    const double a = angle + pi / 4 + i * pi / 2;
    p.push_back({cx + h * std::sqrt(2.0) * std::cos(a), cy + h * std::sqrt(2.0) * std::sin(a), 0});
  }
  if (!ccw) std::reverse(p.begin(), p.end());
  return polygon(p);
}

/// Interval of a horizontal line y = c inside a convex CCW polygon (empty: lo >= hi).
std::pair<double, double> convex_chord(const std::vector<Vec3>& p, double c) {
  double lo = 1e300, hi = -1e300;
  for (std::size_t i = 0; i < p.size(); ++i) {
    const Vec3 a = p[i], b = p[(i + 1) % p.size()];
    if ((a.y - c) * (b.y - c) > 0 || a.y == b.y) continue;
    const double x = a.x + (c - a.y) / (b.y - a.y) * (b.x - a.x);
    lo = std::min(lo, x);
    hi = std::max(hi, x);
  }
  return {lo, hi};
}

double total_length(const Hatching& h) {
  double s = 0;
  for (const auto& g : h.segments) s += g.length();
  return s;
}

double signed_area(const std::vector<NurbsCurve>& cs) {  // Green, Gauss-Legendre per span
  static const double gx[5] = {-0.9061798459386640, -0.5384693101056831, 0.0, 0.5384693101056831,
                               0.9061798459386640};
  static const double gw[5] = {0.2369268850561891, 0.4786286704993665, 0.5688888888888889,
                               0.4786286704993665, 0.2369268850561891};
  double a = 0;
  for (const auto& c : cs) {
    const double u0 = c.domain_start(), u1 = c.domain_end();
    for (int k = 0; k < 400; ++k) {
      const double a0 = u0 + (u1 - u0) * k / 400, a1 = u0 + (u1 - u0) * (k + 1) / 400;
      for (int g = 0; g < 5; ++g) {
        const auto d = c.derivatives(0.5 * (a0 + a1) + 0.5 * (a1 - a0) * gx[g], 1);
        a += 0.25 * (a1 - a0) * gw[g] * (d[0].x * d[1].y - d[0].y * d[1].x);
      }
    }
  }
  return a;
}

/// Contour points lie at distance d from the nearest loop, consecutive curves meet and close.
void check_contour(const std::vector<Loop>& loops, const HatchContour& hc, double d) {
  REQUIRE(!hc.curves.empty());
  for (std::size_t i = 0; i < hc.curves.size(); ++i) {
    const auto& c = hc.curves[i];
    const auto& n = hc.curves[(i + 1) % hc.curves.size()];
    CHECK(distance(c.evaluate(c.domain_end()), n.evaluate(n.domain_start())) < 1e-4);
    for (const auto& p : c.sample(15)) {
      double dist = 1e300;
      for (const auto& l : loops)
        for (const auto& lc : l) dist = std::min(dist, closest_point(lc, p).distance);
      CHECK(dist == Approx(d).margin(1e-6));
    }
  }
}

}  // namespace

TEST_CASE("hatch: chords of a circle are exact", "[hatch]") {
  const double R = 3.0, cx = 0.37, cy = -0.21, h = 0.25;
  const auto res = hatch({circle(cx, cy, R)}, {1, 0, 0}, h);
  // lines y = k h with |y - cy| < R
  const auto k0 = static_cast<std::int64_t>(std::ceil((cy - R) / h)),
             k1 = static_cast<std::int64_t>(std::floor((cy + R) / h));
  REQUIRE(res.segments.size() == static_cast<std::size_t>(k1 - k0 + 1));
  for (const auto& s : res.segments) {
    const double y = static_cast<double>(s.line) * h;
    const double half = std::sqrt(R * R - (y - cy) * (y - cy));
    CHECK(s.start.y == Approx(y).margin(1e-14));
    CHECK(s.end.y == Approx(y).margin(1e-14));
    CHECK(s.start.x == Approx(cx - half).margin(1e-12));
    CHECK(s.end.x == Approx(cx + half).margin(1e-12));
    CHECK(s.region == 0);
  }
  // sorted by line
  for (std::size_t i = 1; i < res.segments.size(); ++i)
    CHECK(res.segments[i].line > res.segments[i - 1].line);
}

TEST_CASE("hatch: rotated square, count and lengths", "[hatch]") {
  const double ang = 0.3;
  std::vector<Vec3> p;
  for (int i = 0; i < 4; ++i) {
    const double a = ang + pi / 4 + i * pi / 2;
    p.push_back({2 * std::sqrt(2.0) * std::cos(a), 2 * std::sqrt(2.0) * std::sin(a), 0});
  }
  const double h = 0.1;
  const auto res = hatch({polygon(p)}, {1, 0, 0}, h, {.phase = 0.05});
  std::size_t expected = 0;
  for (int k = -100; k <= 100; ++k) {
    const double c = 0.05 + k * h;
    const auto [lo, hi] = convex_chord(p, c);
    if (lo >= hi) continue;
    ++expected;
    bool found = false;
    for (const auto& s : res.segments)
      if (s.line == k) {
        found = true;
        CHECK(s.start.x == Approx(lo).margin(1e-12));
        CHECK(s.end.x == Approx(hi).margin(1e-12));
      }
    CHECK(found);
  }
  CHECK(res.segments.size() == expected);
}

TEST_CASE("hatch: spacing times total length approximates the area", "[hatch]") {
  // ring with an off-centre hole, oblique direction
  const double R = 2.0, r = 0.7;
  const double A = pi * (R * R - r * r);
  for (const double h : {0.1, 0.02}) {
    const auto res =
        hatch({circle(0, 0, R), circle(0.4, 0.3, r, false)}, {1, 0.4, 0}, h, {.phase = 0.013});
    const double err = std::abs(h * total_length(res) - A);
    CHECK(err < 2 * pi * (R + r) * h);  // O(h)
  }
}

TEST_CASE("hatch: ring and islands", "[hatch]") {
  SECTION("ring: lines through the hole give two segments") {
    const auto res =
        hatch({circle(0, 0, 2), circle(0, 0, 1, false)}, {1, 0, 0}, 0.5, {.phase = 0.25});
    for (const auto& s : res.segments) {
      const double y = s.start.y;
      CHECK(s.region == 0);
      if (std::abs(y) < 1) {
        // [-sqrt(4-y^2), -sqrt(1-y^2)] and its mirror image
        const double outer = std::max(std::abs(s.start.x), std::abs(s.end.x)),
                     inner = std::min(std::abs(s.start.x), std::abs(s.end.x));
        CHECK(outer == Approx(std::sqrt(4 - y * y)).margin(1e-12));
        CHECK(inner == Approx(std::sqrt(1 - y * y)).margin(1e-12));
        CHECK((s.start_loop == 1) != (s.end_loop == 1));
      }
    }
    std::size_t through_hole = 0;
    for (const auto& s : res.segments) through_hole += std::abs(s.start.y) < 1 ? 1 : 0;
    CHECK(through_hole == 8);  // y = +-0.25, +-0.75: two segments each
  }
  SECTION("two islands and an island in a hole") {
    const auto res =
        hatch({square(-3, 0, 1), circle(3, 0, 2), circle(3, 0, 1, false), circle(3, 0, 0.5)},
              {1, 0, 0}, 0.1, {.phase = 0.05});
    bool left = false, island = false;
    for (const auto& s : res.segments) {
      const double xm = 0.5 * (s.start.x + s.end.x);
      if (xm < 0) {
        CHECK(s.region == 0);
        CHECK(s.length() == Approx(2.0).margin(1e-12));
        left = true;
      } else if (std::abs(xm - 3) < 0.5 && std::abs(s.start.y) < 0.5) {
        CHECK(s.region == 3);  // island inside the hole
        island = true;
      } else {
        CHECK(s.region == 1);
      }
    }
    CHECK(left);
    CHECK(island);
  }
}

TEST_CASE("hatch: tangent lines and edges on lines give no degenerate segments", "[hatch]") {
  SECTION("circle touched by the lines y = +-1") {
    const auto res = hatch({circle(0, 0, 1)}, {1, 0, 0}, 0.5);
    CHECK(res.segments.size() == 3);  // y = -0.5, 0, 0.5
    for (const auto& s : res.segments) CHECK(s.length() > 0.5);
  }
  SECTION("hole touched from outside: the line is not split") {
    const auto res = hatch({circle(0, 0, 3), circle(0, 0, 1, false)}, {1, 0, 0}, 1.0);
    for (const auto& s : res.segments)
      if (std::abs(std::abs(s.start.y) - 1) < 1e-12)
        CHECK(s.length() == Approx(2 * std::sqrt(8.0)).margin(1e-12));
  }
  SECTION("axis-aligned square with edges on hatch lines") {
    const auto res =
        hatch({polygon({{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}})}, {1, 0, 0}, 0.25);
    REQUIRE(res.segments.size() == 3);  // y = 0.25, 0.5, 0.75; y = 0 and y = 1 are edges
    for (const auto& s : res.segments) CHECK(s.length() == Approx(1.0).margin(1e-14));
  }
  SECTION("diamond with vertices on hatch lines") {
    const auto res = hatch({square(0, 0, 1, pi / 4)}, {1, 0, 0}, std::sqrt(2.0) / 4);
    for (const auto& s : res.segments) CHECK(s.length() > 0.1);
    CHECK(res.segments.size() == 7);  // the top and bottom vertices only touch
  }
  SECTION("L-shape: an edge on a line that is crossed") {
    const auto res =
        hatch({polygon({{0, 0, 0}, {2, 0, 0}, {2, 1, 0}, {1, 1, 0}, {1, 2, 0}, {0, 2, 0}})},
              {1, 0, 0}, 0.5);
    for (const auto& s : res.segments)
      if (s.line == 2) {  // y = 1: the section from above, x in [0, 1]
        CHECK(s.start.x == Approx(0.0).margin(1e-14));
        CHECK(s.end.x == Approx(1.0).margin(1e-14));
      }
  }
}

TEST_CASE("hatch: lines satisfy n . p = phase + k spacing for any direction", "[hatch]") {
  const Vec3 d = normalized(Vec3{0.6, -0.8, 0});
  const Vec3 n{-d.y, d.x, 0};
  const double h = 0.3, ph = 0.11;
  const auto res = hatch({circle(1, 2, 2.5)}, d, h, {.phase = ph, .min_length = 0.5});
  REQUIRE(!res.segments.empty());
  for (const auto& s : res.segments) {
    const double c = ph + static_cast<double>(s.line) * h;
    CHECK(dot(n, s.start) == Approx(c).margin(1e-12));
    CHECK(dot(n, s.end) == Approx(c).margin(1e-12));
    CHECK(dot(d, s.end - s.start) > 0.5);  // along d, min_length respected
    CHECK(distance(s.start, Vec3{1, 2, 0}) == Approx(2.5).margin(1e-12));
  }
}

TEST_CASE("hatch: contours", "[hatch]") {
  SECTION("circle: concentric circle") {
    const auto res = hatch({circle(0.5, 0, 2)}, {1, 0, 0}, 0.5, {.contours = {0.25, 0.5}});
    REQUIRE(res.contours.size() == 2);
    CHECK(res.contours[1].distance == 0.5);
    for (const auto& hc : res.contours)
      for (const auto& c : hc.curves)
        for (const auto& p : c.sample(50))
          CHECK(distance(p, Vec3{0.5, 0, 0}) == Approx(2 - hc.distance).margin(1e-6));
  }
  SECTION("square: trimmed at convex corners") {
    const auto res = hatch({square(0, 0, 1)}, {1, 0, 0}, 0.5, {.contours = {0.25}});
    REQUIRE(res.contours.size() == 1);
    const auto& cs = res.contours[0].curves;
    REQUIRE(cs.size() == 4);
    for (const auto& c : cs) {
      const Vec3 a = c.evaluate(c.domain_start()), b = c.evaluate(c.domain_end());
      CHECK(std::abs(a.x) == Approx(0.75).margin(1e-6));
      CHECK(std::abs(a.y) == Approx(0.75).margin(1e-6));
      CHECK(distance(a, b) == Approx(1.5).margin(1e-6));
    }
  }
  SECTION("square hole: arcs at the (reflex) corners") {
    const auto res =
        hatch({square(0, 0, 3), square(0, 0, 1, 0.0, false)}, {1, 0, 0}, 0.5, {.contours = {0.25}});
    REQUIRE(res.contours.size() == 2);
    const auto& hole = res.contours[1];
    CHECK(hole.loop == 1);
    REQUIRE(hole.curves.size() == 8);  // 4 sides + 4 arcs
    double len = 0;
    for (const auto& c : hole.curves) len += arc_length(c);
    CHECK(len == Approx(8.0 + 2 * pi * 0.25).epsilon(1e-6));
    // the contour runs at distance 0.25 from the square hole, outside it
    for (const auto& c : hole.curves)
      for (const auto& p : c.sample(20)) {
        const double dx = std::max(std::abs(p.x) - 1, 0.0), dy = std::max(std::abs(p.y) - 1, 0.0);
        CHECK(std::hypot(dx, dy) == Approx(0.25).margin(1e-6));
      }
  }
  SECTION("vanishing contours") {
    CHECK(hatch({square(0, 0, 1)}, {1, 0, 0}, 0.5, {.contours = {1.2}}).contours.empty());
    CHECK(hatch({circle(0, 0, 1)}, {1, 0, 0}, 0.5, {.contours = {1.5}}).contours.empty());
    // island survives at 0.2, vanishes at 0.6; the outer circle keeps its contour
    const auto res =
        hatch({circle(0, 0, 3), square(0, 0, 0.5)}, {1, 0, 0}, 0.5, {.contours = {0.2, 0.6}});
    REQUIRE(res.contours.size() == 3);
    CHECK(res.contours[2].loops == std::vector<std::size_t>{1});
    CHECK(res.contours[2].index == 0);
  }
}

TEST_CASE("hatch: contours split, merge and lose cusps", "[hatch]") {
  SECTION("two holes 0.4 apart merge at 0.3") {
    const std::vector<Loop> loops{circle(0, 0, 5), circle(-1.2, 0, 1, false),
                                  circle(1.2, 0, 1, false)};
    const auto res = hatch(loops, {1, 0, 0}, 0.5, {.contours = {0.3}});
    REQUIRE(res.contours.size() == 2);
    CHECK(res.contours[0].loops == std::vector<std::size_t>{0});
    CHECK(res.contours[1].loops == std::vector<std::size_t>{1, 2});
    for (const auto& hc : res.contours) check_contour(loops, hc, 0.3);
    // the merged hole contour: two arcs of radius 1.3 meeting on the axis (area by Green)
    const double a = std::acos(1.2 / 1.3);  // half angle of the arcs at the waist
    const double area = 2 * ((pi - a) * 1.3 * 1.3 + 1.2 * 1.3 * std::sin(a));
    CHECK(std::abs(signed_area(res.contours[1].curves)) == Approx(area).epsilon(1e-7));
  }
  SECTION("hole close to the outer loop merges with it") {
    const std::vector<Loop> loops{circle(0, 0, 5), circle(3.7, 0, 1, false)};
    const auto res = hatch(loops, {1, 0, 0}, 0.5, {.contours = {0.25}});
    REQUIRE(res.contours.size() == 1);
    CHECK(res.contours[0].loops == std::vector<std::size_t>{0, 1});
    check_contour(loops, res.contours[0], 0.25);
  }
  SECTION("dumbbell: a neck narrower than 2 d splits the contour") {
    const std::vector<Loop> loops{polygon({{-3, -1, 0},
                                           {-1, -1, 0},
                                           {-1, -0.2, 0},
                                           {1, -0.2, 0},
                                           {1, -1, 0},
                                           {3, -1, 0},
                                           {3, 1, 0},
                                           {1, 1, 0},
                                           {1, 0.2, 0},
                                           {-1, 0.2, 0},
                                           {-1, 1, 0},
                                           {-3, 1, 0}})};
    const auto one = hatch(loops, {1, 0, 0}, 0.5, {.contours = {0.1}});
    REQUIRE(one.contours.size() == 1);
    check_contour(loops, one.contours[0], 0.1);
    const auto two = hatch(loops, {1, 0, 0}, 0.5, {.contours = {0.3}});
    REQUIRE(two.contours.size() == 2);
    for (const auto& hc : two.contours) {
      CHECK(hc.loops == std::vector<std::size_t>{0});
      // the lobe shrunk to a square of side 1.4, bulging into the neck between arcs of radius 0.3
      // about the reflex corners: 3 + 2 sides, 2 arcs
      CHECK(hc.curves.size() == 7);
      const double bulge = 2 * (0.06 - 0.5 * (0.2 * std::sqrt(0.05) + 0.09 * std::asin(2.0 / 3.0)));
      CHECK(signed_area(hc.curves) == Approx(1.4 * 1.4 + bulge).epsilon(1e-8));
      check_contour(loops, hc, 0.3);
    }
  }
  SECTION("ellipse: the offset beyond the radius of curvature loses its swallowtails") {
    // a = 3, b = 1: radius of curvature 1/3 at the tips, d = 0.5
    const auto e = make_circle({0, 0, 0}, {0, 0, 1}, 1)
                       .transformed(Transform::scaling({0, 0, 0}, Vec3{3, 1, 1}));
    const std::vector<Loop> loops{{e}};
    const auto res = hatch(loops, {1, 0, 0}, 0.5, {.contours = {0.5}});
    REQUIRE(res.contours.size() == 1);
    check_contour(loops, res.contours[0], 0.5);
    CHECK(res.contours[0].curves.size() >= 2);  // trimmed where the two halves cross
  }
}

TEST_CASE("hatch: section loops of a shell", "[hatch]") {
  // cylinder r = 10, z = 0..5, cut at z = 2.5: one exact circle in the plane coordinates
  const double R = 10;
  const auto side = extrude(make_circle({0, 0, 0}, {0, 0, 1}, R), {0, 0, 5});
  const auto lid = [&](double z) {
    return TrimmedSurface(make_rectangle({-R, -R, z}, {2 * R, 0, 0}, {0, 2 * R, 0}),
                          circle_loop(0.5, 0.5, 0.5));
  };
  const auto shell = Shell::sew({TrimmedSurface(side), lid(0), lid(5)});
  const auto sec = section(shell, Plane{{0, 0, 2.5}, {0, 0, 1}});
  REQUIRE(sec.loops.size() == 1);
  const Vec3 axis{0, 0, 2.5};
  const Vec3 c{dot(axis - sec.origin, sec.s_axis), dot(axis - sec.origin, sec.t_axis), 0};
  const auto res = hatch(sec.loops, {1, 0, 0}, 1.0, {.phase = 0.5, .contours = {0.5}});
  REQUIRE(!res.segments.empty());
  for (const auto& s : res.segments) {
    const double y = s.start.y - c.y;
    CHECK(s.length() == Approx(2 * std::sqrt(R * R - y * y)).margin(1e-6));
  }
  REQUIRE(res.contours.size() == 1);
  for (const auto& cv : res.contours[0].curves)
    for (const auto& p : cv.sample(20)) CHECK(distance(p, c) == Approx(R - 0.5).margin(1e-5));
}

TEST_CASE("hatch: argument checks", "[hatch]") {
  CHECK_THROWS_AS(hatch({circle(0, 0, 1)}, {1, 0, 0}, 0.0), GeometryError);
  CHECK_THROWS_AS(hatch({circle(0, 0, 1)}, {0, 0, 1}, 0.1), GeometryError);
  CHECK_THROWS_AS(hatch({circle(0, 0, 1)}, {1, 0, 0}, 1e-9), GeometryError);
  CHECK_THROWS_AS(hatch({circle(0, 0, 1)}, {1, 0, 0}, 0.1, {.contours = {-0.1}}), GeometryError);
  CHECK(hatch(std::vector<std::vector<NurbsCurve>>{}, {1, 0, 0}, 0.1).segments.empty());
}
