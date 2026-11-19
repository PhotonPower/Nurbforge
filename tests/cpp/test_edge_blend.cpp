#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <numbers>

#include "nurbforge/core/error.hpp"
#include "nurbforge/topology/boolean.hpp"
#include "nurbforge/topology/edge_blend.hpp"
#include "nurbforge/topology/solids.hpp"

using namespace nf;
constexpr double pi = std::numbers::pi;

namespace {

// Volume of the ring swept by the region between a right-angled edge and a fillet of radius r
// (square r x r minus the quarter disc) at distance R from the axis (Pappus): `inward` for an
// edge on the outside of the ring (the region lies at R - x), outward otherwise.
double fillet_ring(double R, double r, bool inward) {
  const double A = r * r * (1 - pi / 4);
  const double xq = r - 4 * r / (3 * pi);  // centroid of the quarter disc from the corner
  const double xbar = (r * r * r / 2 - pi * r * r / 4 * xq) / A;
  return 2 * pi * (inward ? R - xbar : R + xbar) * A;
}

std::size_t count_label(const Shell& s, const std::string& l) { return s.faces_labeled(l).size(); }

}  // namespace

TEST_CASE("fillet and chamfer of box edges", "[edge_blend]") {
  const Shell box = make_box({0, 0, 0}, {4, 3, 2});
  const auto e = edges_between(box, "top", "front");
  REQUIRE(e.size() == 1);
  const double r = 0.5;
  SECTION("fillet: exact volume, faces trimmed, end faces rerouted") {
    const Shell f = fillet_edges(box, e, r);
    CHECK(f.is_closed());
    CHECK(f.is_consistently_oriented());
    CHECK(f.faces().size() == 7);
    CHECK(count_label(f, "fillet") == 1);
    CHECK(std::abs(f.volume() - (24 - 4 * r * r * (1 - pi / 4))) <= 1e-7);
    CHECK(std::abs(f.area() - (2 * (12 + 8 + 6) - 2 * 4 * r - 2 * r * r * (1 - pi / 4) +
                               4 * pi * r / 2)) <= 1e-6);
  }
  SECTION("chamfer and custom label") {
    const Shell c = chamfer_edges(box, e, r, {"bevel", 0.0});
    CHECK(c.is_closed());
    CHECK(count_label(c, "bevel") == 1);
    CHECK(std::abs(c.volume() - (24 - 4 * r * r / 2)) <= 1e-9);
  }
  SECTION("two edges without a common vertex") {
    auto two = edges_between(box, "top", "front");
    for (const auto i : edges_between(box, "bottom", "back")) two.push_back(i);
    const Shell f = fillet_edges(box, two, r);
    CHECK(std::abs(f.volume() - (24 - 8 * r * r * (1 - pi / 4))) <= 1e-7);
  }
  SECTION("errors") {
    auto corner = edges_between(box, "top", "front");
    for (const auto i : edges_between(box, "top", "left")) corner.push_back(i);
    CHECK_THROWS_AS(fillet_edges(box, corner, r), GeometryError);  // vertex blend
    CHECK_THROWS_AS(fillet_edges(box, e, 0.0), GeometryError);
    CHECK_THROWS_AS(fillet_edges(box, e, 2.5), GeometryError);  // larger than the faces
    CHECK_THROWS_AS(chamfer_edges(box, {99}, 0.1), GeometryError);
    CHECK(fillet_edges(box, {}, r).faces().size() == 6);
  }
}

TEST_CASE("fillet and chamfer of closed edges", "[edge_blend]") {
  const double R = 2.0, r = 0.4;
  const Shell cyl = make_cylinder({0, 0, 0}, {0, 0, 1}, R, 3.0);
  const auto rim = edges_between(cyl, "top", "side");
  SECTION("rim of a cylinder") {
    const Shell f = fillet_edges(cyl, rim, r);
    CHECK(f.is_closed());
    CHECK(std::abs(f.volume() - (pi * R * R * 3 - fillet_ring(R, r, true))) <= 1e-6);
    const Shell c = chamfer_edges(cyl, rim, r);
    CHECK(std::abs(c.volume() - (pi * R * R * 3 - 2 * pi * (R - r / 3) * r * r / 2)) <= 1e-6);
  }
  SECTION("bore through a block (convex) and base of a boss (concave: adds material)") {
    const Shell block =
        boolean(make_box({-3, -3, 0}, {6, 6, 2}), make_cylinder({0, 0, -1}, {0, 0, 1}, 1.0, 4.0),
                BooleanOp::Difference);
    const Shell f = fillet_edges(block, edges_between(block, "top", "side"), 0.3);
    CHECK(std::abs(f.volume() - (block.volume() - fillet_ring(1.0, 0.3, false))) <= 1e-6);
    const Shell boss = boolean(make_box({-3, -3, -2}, {6, 6, 2}),
                               make_cylinder({0, 0, -1}, {0, 0, 1}, 1.0, 3.0), BooleanOp::Union);
    // the base circle (box top / cylinder side); the cylinder's own top rim is excluded
    std::vector<std::size_t> base;
    for (const auto i : edges_between(boss, "top", "side"))
      if (std::abs(boss.edges()[i].curve.evaluate(boss.edges()[i].curve.domain_start()).z) < 1e-6)
        base.push_back(i);
    REQUIRE(!base.empty());
    const Shell c = chamfer_edges(boss, base, 0.3);
    CHECK(c.is_closed());
    CHECK(std::abs(c.volume() - (boss.volume() + 2 * pi * (1.0 + 0.1) * 0.045)) <= 1e-6);
    const Shell g = fillet_edges(boss, base, 0.3);
    CHECK(std::abs(g.volume() - (boss.volume() + fillet_ring(1.0, 0.3, false))) <= 1e-6);
  }
}
