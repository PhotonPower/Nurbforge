#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <numbers>

#include "nurbforge/geometry/primitives.hpp"
#include "nurbforge/topology/region.hpp"

using namespace nf;
constexpr double pi = std::numbers::pi;

namespace {

bool close(double a, double b, double rel = 1e-9) { return std::abs(a - b) <= rel * std::abs(b); }
Region square(double x, double y, double s) {
  return {{make_rectangle_curve({x, y, 0}, {s, 0, 0}, {0, s, 0})}};
}
Region disc(double x, double y, double r) { return {{make_circle({x, y, 0}, {0, 0, 1}, r)}}; }

}  // namespace

TEST_CASE("region Booleans of squares and discs", "[region]") {
  SECTION("overlapping squares") {
    const Region a = square(0, 0, 2), b = square(1, 1, 2);
    CHECK(close(region_area(region_boolean(a, b, BooleanOp::Union)), 7.0));
    CHECK(close(region_area(region_boolean(a, b, BooleanOp::Intersection)), 1.0));
    CHECK(close(region_area(region_boolean(a, b, BooleanOp::Difference)), 3.0));
    CHECK(region_boolean(a, b, BooleanOp::Union).size() == 1);
  }
  SECTION("disc minus a smaller disc: a ring (outer loop and hole)") {
    const Region ring = region_boolean(disc(0, 0, 2), disc(0, 0, 1), BooleanOp::Difference);
    REQUIRE(ring.size() == 2);
    CHECK(close(region_area(ring), 3 * pi, 1e-9));
    // the hole is clockwise
    double neg = 0;
    for (const auto& l : ring) neg += TrimLoop{l}.signed_area() < 0 ? 1 : 0;
    CHECK(neg == 1);
  }
  SECTION("lens of two discs, and a disc cut out of a square edge") {
    const double d = 1.0;  // two unit discs at distance 1: lens area 2 acos(d/2) - d/2 sqrt(4-d^2)
    const double lens = 2 * std::acos(d / 2) - d / 2 * std::sqrt(4 - d * d);
    CHECK(close(region_area(region_boolean(disc(0, 0, 1), disc(d, 0, 1), BooleanOp::Intersection)),
                lens, 1e-7));
    CHECK(close(region_area(region_boolean(disc(0, 0, 1), disc(d, 0, 1), BooleanOp::Union)),
                2 * pi - lens, 1e-7));
    CHECK(
        close(region_area(region_boolean(square(0, 0, 2), disc(0, 1, 0.5), BooleanOp::Difference)),
              4 - pi * 0.25 / 2, 1e-7));
  }
  SECTION("touching and disjoint squares") {
    const Region u = region_boolean(square(0, 0, 1), square(1, 0, 1), BooleanOp::Union);
    CHECK(u.size() == 1);  // the shared edge cancels
    CHECK(close(region_area(u), 2.0));
    const Region two = region_boolean(square(0, 0, 1), square(3, 0, 1), BooleanOp::Union);
    CHECK(two.size() == 2);
    CHECK(close(region_area(two), 2.0));
    CHECK(region_boolean(square(0, 0, 1), square(3, 0, 1), BooleanOp::Intersection).empty());
  }
  SECTION("input orientation does not matter; a square with a hole") {
    Region holed = square(0, 0, 4);
    holed.push_back({make_circle({2, 2, 0}, {0, 0, 1}, 1.0)});  // same orientation as the outer
    CHECK(close(region_area(holed), 16 - pi, 1e-12));
    const Region filled = region_boolean(holed, disc(2, 2, 1.5), BooleanOp::Union);
    CHECK(filled.size() == 1);
    CHECK(close(region_area(filled), 16.0, 1e-9));
  }
}
