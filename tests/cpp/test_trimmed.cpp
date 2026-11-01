#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <numbers>
#include <vector>

#include "nurbforge/analysis/properties.hpp"
#include "nurbforge/core/error.hpp"
#include "nurbforge/geometry/primitives.hpp"
#include "nurbforge/modeling/construct.hpp"
#include "nurbforge/topology/trimmed_surface.hpp"

using namespace nf;
constexpr double pi = std::numbers::pi;

namespace {

NurbsSurface plate() { return make_rectangle({0, 0, 0}, {4, 0, 0}, {0, 2, 0}); }  // x = 4u, y = 2v

NurbsSurface sphere(double R) {
  return revolve(make_arc({0, 0, 0}, {0, 0, -1}, {1, 0, 0}, R, 0.0, pi), {0, 0, 0}, {0, 0, 1}, 2 * pi);
}

bool rel_close(double a, double b, double rel) { return std::abs(a - b) <= rel * std::abs(b); }

}  // namespace

TEST_CASE("untrimmed surface", "[trim]") {
  const TrimmedSurface t(plate());
  CHECK(t.contains(0.5, 0.5));
  CHECK_FALSE(t.contains(1.5, 0.5));
  CHECK(rel_close(t.area(), 8.0, 1e-10));
  CHECK(t.holes().empty());
}

TEST_CASE("circular aperture, hole and nested loops on a plate", "[trim]") {
  const auto s = plate();
  const TrimmedSurface disk(s, circle_loop(0.5, 0.5, 0.25));  // ellipse with semi-axes 1 and 0.5
  CHECK(disk.contains(0.5, 0.5));
  CHECK(disk.contains(0.7, 0.5));
  CHECK_FALSE(disk.contains(0.8, 0.5));
  CHECK(rel_close(disk.area(), pi * 1.0 * 0.5, 1e-5));

  const auto holed = TrimmedSurface(s).with_hole(circle_loop(0.5, 0.5, 0.25));
  CHECK_FALSE(holed.contains(0.5, 0.5));
  CHECK(holed.contains(0.1, 0.1));
  CHECK(rel_close(holed.area(), 8.0 - pi * 0.5, 1e-5));

  // the 3D boundary of the disk is the ellipse (x-2)^2 + 4 (y-1)^2 = 1
  const auto b = disk.boundary_curves(1e-8);
  REQUIRE(b.size() == 1);
  for (const auto& c : b[0])
    for (const auto& p : c.sample(100)) CHECK(std::abs((p.x - 2) * (p.x - 2) + 4 * (p.y - 1) * (p.y - 1) - 1) < 1e-7);
}

TEST_CASE("loop validation and orientation", "[trim]") {
  const auto s = plate();
  // clockwise outer loop is reversed on construction
  const TrimmedSurface t(s, circle_loop(0.5, 0.5, 0.25).reversed());
  CHECK(t.outer().signed_area() > 0.0);
  const auto h = TrimmedSurface(s).with_hole(circle_loop(0.5, 0.5, 0.1));
  CHECK(h.holes()[0].signed_area() < 0.0);
  TrimLoop open{{make_line({0.1, 0.1, 0}, {0.9, 0.1, 0}), make_line({0.9, 0.1, 0}, {0.9, 0.9, 0})}};
  CHECK_THROWS_AS(TrimmedSurface(s, open), GeometryError);
  CHECK_THROWS_AS(TrimmedSurface(s, circle_loop(0.9, 0.5, 0.25)), GeometryError);  // leaves the domain
  CHECK(rel_close(TrimmedSurface(s, polygon_loop({{0.1, 0.1}, {0.6, 0.1}, {0.1, 0.6}})).area(), 0.5 * 2.0 * 1.0, 1e-5));
}

TEST_CASE("split along chords and nested closed curves", "[trim]") {
  const auto s = make_rectangle({0, 0, 0}, {1, 0, 0}, {0, 1, 0});  // x = u, y = v
  // two non-crossing chords: x = 0.3 and x = 0.7 -> three strips
  const std::vector<NurbsCurve> chords{make_line({0.3, 0, 0}, {0.3, 1, 0}), make_line({0.7, 1, 0}, {0.7, 0, 0})};
  auto parts = split(s, chords);
  REQUIRE(parts.size() == 3);
  std::vector<double> areas;
  for (const auto& p : parts) areas.push_back(p.area());
  std::sort(areas.begin(), areas.end());
  CHECK(rel_close(areas[0], 0.3, 1e-8));
  CHECK(rel_close(areas[1], 0.3, 1e-8));
  CHECK(rel_close(areas[2], 0.4, 1e-8));

  // two concentric circles: outside, annulus, disk
  const std::vector<NurbsCurve> rings{make_circle({0.5, 0.5, 0}, {0, 0, 1}, 0.4), make_circle({0.5, 0.5, 0}, {0, 0, 1}, 0.2)};
  parts = split(s, rings);
  REQUIRE(parts.size() == 3);
  areas.clear();
  for (const auto& p : parts) areas.push_back(p.area());
  std::sort(areas.begin(), areas.end());
  CHECK(rel_close(areas[0], pi * 0.04, 1e-5));
  CHECK(rel_close(areas[1], pi * (0.16 - 0.04), 1e-5));
  CHECK(rel_close(areas[2], 1.0 - pi * 0.16, 1e-5));
  CHECK_THROWS_AS(split(s, std::vector<NurbsCurve>{make_line({0.2, 0.2, 0}, {0.5, 0.5, 0})}), GeometryError);
}

TEST_CASE("split a plate by a sphere and a sphere by a plane", "[trim][ssi]") {
  const double R = 2.0, h = 0.8;
  const auto big = make_rectangle({-3, -3, h}, {6, 0, 0}, {0, 6, 0});
  const auto ball = sphere(R);
  const auto plate_parts = split(big, ball);
  REQUIRE(plate_parts.size() == 2);
  const double rho2 = R * R - h * h;
  std::vector<double> a{plate_parts[0].area(1e-7), plate_parts[1].area(1e-7)};
  std::sort(a.begin(), a.end());
  CHECK(rel_close(a[0], pi * rho2, 1e-5));
  CHECK(rel_close(a[1], 36.0 - pi * rho2, 1e-5));

  // the circle on the sphere crosses its seam: open chord from u = 0 to u = 1 -> two caps
  const auto caps = split(ball, big);
  REQUIRE(caps.size() == 2);
  std::vector<double> c{caps[0].area(1e-7), caps[1].area(1e-7)};
  std::sort(c.begin(), c.end());
  CHECK(rel_close(c[0], 2 * pi * R * (R - h), 1e-5));
  CHECK(rel_close(c[1], 2 * pi * R * (R + h), 1e-5));
}

TEST_CASE("general split of trimmed faces (planar arrangement)", "[trim][arrangement]") {
  const auto sq = make_rectangle({0, 0, 0}, {1, 0, 0}, {0, 1, 0});  // x = u, y = v
  const auto areas = [](const std::vector<TrimmedSurface>& parts) {
    std::vector<double> a;
    for (const auto& p : parts) a.push_back(p.area(1e-9));
    std::sort(a.begin(), a.end());
    return a;
  };
  // crossing chords -> four quarters
  auto a = areas(split(TrimmedSurface(sq), {make_line({0.5, -0.2, 0}, {0.5, 1.2, 0}), make_line({-0.2, 0.5, 0}, {1.2, 0.5, 0})}));
  REQUIRE(a.size() == 4);
  for (double x : a) CHECK(std::abs(x - 0.25) < 1e-9);

  // disk with a chord sticking out of it: two circular segments
  const TrimmedSurface disk(sq, circle_loop(0.5, 0.5, 0.4));
  a = areas(split(disk, {make_line({0.7, -1, 0}, {0.7, 2, 0})}));
  REQUIRE(a.size() == 2);
  const double r = 0.4, d = 0.2, seg = r * r * std::acos(d / r) - d * std::sqrt(r * r - d * d);
  CHECK(rel_close(a[0], seg, 1e-6));
  CHECK(rel_close(a[1], pi * r * r - seg, 1e-6));

  // square with a hole, cut through the hole: two halves with a notch each
  const auto holed = TrimmedSurface(sq).with_hole(circle_loop(0.5, 0.5, 0.2));
  a = areas(split(holed, {make_line({0.5, -0.5, 0}, {0.5, 1.5, 0})}));
  REQUIRE(a.size() == 2);
  for (double x : a) CHECK(rel_close(x, 0.5 * (1 - pi * 0.04), 1e-6));

  // closed cut curve: inside and outside
  a = areas(split(TrimmedSurface(sq), {make_circle({0.4, 0.6, 0}, {0, 0, 1}, 0.25)}));
  REQUIRE(a.size() == 2);
  CHECK(rel_close(a[0], pi * 0.0625, 1e-6));
  CHECK(rel_close(a[1], 1 - pi * 0.0625, 1e-6));

  // T-junction: a vertical cut from the top ends on a horizontal cut -> three regions
  a = areas(split(TrimmedSurface(sq), {make_line({0, 0.4, 0}, {1, 0.4, 0}), make_line({0.3, 1, 0}, {0.3, 0.4, 0})}));
  REQUIRE(a.size() == 3);
  CHECK(std::abs(a[0] - 0.18) < 1e-9);
  CHECK(std::abs(a[1] - 0.4) < 1e-9);
  CHECK(std::abs(a[2] - 0.42) < 1e-9);
  // a dangling cut (ends inside) changes nothing
  CHECK(split(TrimmedSurface(sq), {make_line({0.2, 0.2, 0}, {0.6, 0.5, 0})}).size() == 1);
}
