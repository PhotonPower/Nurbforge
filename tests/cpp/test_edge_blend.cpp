#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <numbers>
#include <vector>

#include "nurbforge/core/error.hpp"
#include "nurbforge/geometry/primitives.hpp"
#include "nurbforge/modeling/construct.hpp"
#include "nurbforge/topology/boolean.hpp"
#include "nurbforge/topology/edge_blend.hpp"
#include "nurbforge/topology/solid_ops.hpp"
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

// Corner of a box with its three edges rounded (eighth ball) instead of sharp, minus what the
// three edge fillets alone (run to the corner) would remove there.
double rounded_corner(double r) { return r * r * r * (2 - 7 * pi / 12); }

std::vector<std::size_t> edges_of(const Shell& s, std::initializer_list<const char*> pairs) {
  std::vector<std::size_t> out;
  const std::vector<const char*> p(pairs);
  for (std::size_t i = 0; i + 1 < p.size(); i += 2)
    for (const auto e : edges_between(s, p[i], p[i + 1])) out.push_back(e);
  return out;
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

TEST_CASE("vertex blends", "[edge_blend]") {
  const Shell box = make_box({0, 0, 0}, {4, 3, 2});
  const double r = 0.5, A = r * r * (1 - pi / 4);
  SECTION("two edges at a corner: mitre (fillet and chamfer)") {
    const auto e = edges_of(box, {"top", "front", "top", "left"});
    const Shell f = fillet_edges(box, e, r);
    CHECK(f.is_closed());
    CHECK(f.is_consistently_oriented());
    CHECK(std::abs(f.volume() - (24 - 7 * A + r * r * r * (5.0 / 3 - pi / 2))) <= 1e-7);
    const Shell c = chamfer_edges(box, e, r);
    CHECK(std::abs(c.volume() - (24 - 7 * r * r / 2 + r * r * r / 3)) <= 1e-9);
  }
  SECTION("three edges at a corner: spherical / planar corner patch") {
    const auto e = edges_of(box, {"top", "front", "top", "left", "front", "left"});
    const Shell f = fillet_edges(box, e, r);
    CHECK(f.is_closed());
    CHECK(f.is_consistently_oriented());
    CHECK(count_label(f, "fillet") == 4);
    CHECK(std::abs(f.volume() - (24 - 9 * A + rounded_corner(r))) <= 1e-7);
    const Shell c = chamfer_edges(box, e, r);
    CHECK(count_label(c, "chamfer") == 4);
    CHECK(std::abs(c.volume() - (24 - 9 * r * r / 2 + 2 * r * r * r / 3)) <= 1e-9);
  }
  SECTION("all twelve edges of a box") {
    std::vector<std::size_t> all(box.edges().size());
    for (std::size_t i = 0; i < all.size(); ++i) all[i] = i;
    const Shell f = fillet_edges(box, all, r);
    CHECK(f.is_closed());
    CHECK(f.faces().size() == 6 + 12 + 8);
    CHECK(std::abs(f.volume() - (24 - 36 * A + 8 * rounded_corner(r))) <= 1e-6);
    const Shell c = chamfer_edges(box, all, r);
    CHECK(std::abs(c.volume() - (24 - 36 * r * r / 2 + 8 * 2 * r * r * r / 3)) <= 1e-9);
  }
  SECTION("concave: floor and wall edges of a pocket") {
    const double a = 2, b = 1.5, c = 1;  // pocket size
    const Shell block = boolean(make_box({0, 0, 0}, {4, 3, 2}), make_box({1, 0.75, 1}, {a, b, 3}),
                                BooleanOp::Difference);
    // pocket edges: those between faces whose midpoints lie inside the pocket box
    std::vector<std::size_t> floor, walls;
    for (std::size_t i = 0; i < block.edges().size(); ++i) {
      const auto& k = block.edges()[i].curve;
      const Vec3 m = k.evaluate(0.5 * (k.domain_start() + k.domain_end()));
      const bool in = m.x > 0.99 && m.x < 3.01 && m.y > 0.74 && m.y < 2.26;
      if (!in) continue;
      if (std::abs(m.z - 1) < 1e-9)
        floor.push_back(i);
      else if (m.z > 1 && m.z < 2)
        walls.push_back(i);
    }
    REQUIRE(floor.size() == 4);
    REQUIRE(walls.size() == 4);
    const double V0 = block.volume();
    const Shell f = fillet_edges(block, floor, r);
    CHECK(f.is_closed());
    CHECK(std::abs(f.volume() - (V0 + A * 2 * (a + b) - 4 * r * r * r * (5.0 / 3 - pi / 2))) <=
          1e-7);
    auto both = floor;
    both.insert(both.end(), walls.begin(), walls.end());
    const Shell g = fillet_edges(block, both, r);
    CHECK(g.is_closed());
    CHECK(g.is_consistently_oriented());
    CHECK(std::abs(g.volume() - (V0 + A * (2 * (a + b) + 4 * c) - 4 * rounded_corner(r))) <= 1e-7);
  }
  SECTION("tangent chain: rim of a plate with rounded corners") {
    const double S = 4, R = 1, H = 1, rr = 0.3;
    // outline: lines and quarter arcs (exact), extruded to walls and capped
    std::vector<NurbsCurve> outline;
    outline.push_back(make_line({R, 0, 0}, {S - R, 0, 0}));
    outline.push_back(make_arc({S - R, R, 0}, {1, 0, 0}, {0, 1, 0}, R, -pi / 2, 0));
    outline.push_back(make_line({S, R, 0}, {S, S - R, 0}));
    outline.push_back(make_arc({S - R, S - R, 0}, {1, 0, 0}, {0, 1, 0}, R, 0, pi / 2));
    outline.push_back(make_line({S - R, S, 0}, {R, S, 0}));
    outline.push_back(make_arc({R, S - R, 0}, {1, 0, 0}, {0, 1, 0}, R, pi / 2, pi));
    outline.push_back(make_line({0, S - R, 0}, {0, R, 0}));
    outline.push_back(make_arc({R, R, 0}, {1, 0, 0}, {0, 1, 0}, R, pi, 3 * pi / 2));
    std::vector<TrimmedSurface> walls;
    for (const auto& c : outline)
      walls.push_back(TrimmedSurface(extrude(c, {0, 0, H})).with_label("wall"));
    const Shell plate = cap(Shell::sew(walls));
    const double V0 = H * (S * S - (4 - pi) * R * R);
    REQUIRE(std::abs(std::abs(plate.volume()) - V0) <= 1e-9);
    std::vector<std::size_t> rim;  // edges between the top cap and the walls
    for (const auto e : edges_between(plate, "cap", "wall"))
      if (std::abs(plate.edges()[e].curve.evaluate(plate.edges()[e].curve.domain_start()).z - H) <
          1e-9)
        rim.push_back(e);
    REQUIRE(rim.size() == 8);
    const double Ar = rr * rr * (1 - pi / 4);
    const Shell f = fillet_edges(plate, rim, rr);
    CHECK(f.is_closed());
    CHECK(f.is_consistently_oriented());
    CHECK(count_label(f, "fillet") == 8);
    CHECK(std::abs(f.volume() - (V0 - Ar * 4 * (S - 2 * R) - fillet_ring(R, rr, true))) <= 1e-6);
    CHECK_THROWS_AS(fillet_edges(plate, {rim[0]}, rr), GeometryError);  // part of a tangent chain
    const Shell c = chamfer_edges(plate, rim, rr);
    CHECK(std::abs(c.volume() - (V0 - rr * rr / 2 * 4 * (S - 2 * R) -
                                 2 * pi * (R - rr / 3) * rr * rr / 2)) <= 1e-6);
  }
}

TEST_CASE("variable size", "[edge_blend]") {
  const Shell box = make_box({0, 0, 0}, {4, 3, 2});
  const auto e = edges_between(box, "top", "front");
  const double r0 = 0.3, r1 = 0.6, dr = r1 - r0;
  // integral over the edge (length 4) of size(x)^2 with the smooth step from r0 to r1
  const double I = 4 * (r0 * r0 + r0 * dr + dr * dr * 181.0 / 462);
  EdgeBlendOptions o;
  o.sizes = {{{0, 0, 2}, r0}, {{4, 0, 2}, r1}};
  SECTION("fillet: sections are quarter circles of radius r(x)") {
    const Shell f = fillet_edges(box, e, 1.0, o);  // the size argument is not used
    CHECK(f.is_closed());
    CHECK(f.is_consistently_oriented());
    CHECK(std::abs(f.volume() - (24 - (1 - pi / 4) * I)) <= 1e-7);
  }
  SECTION("chamfer") {
    const Shell c = chamfer_edges(box, e, 1.0, o);
    CHECK(std::abs(c.volume() - (24 - I / 2)) <= 1e-7);
  }
  SECTION("closed chain and errors") {
    const Shell cyl = make_cylinder({0, 0, 0}, {0, 0, 1}, 2.0, 3.0);
    EdgeBlendOptions oc;
    oc.sizes = {{{2, 0, 3}, 0.2}, {{-2, 0, 3}, 0.5}};
    const Shell f = fillet_edges(cyl, edges_between(cyl, "top", "side"), 1.0, oc);
    CHECK(f.is_closed());
    // sections in the planes through the axis: Pappus per angle, r(phi) the smooth step from
    // 0.2 at phi = 0 to 0.5 at phi = pi and back (Simpson)
    const auto rr = [](double phi) {
      double t = phi <= pi ? phi / pi : (2 * pi - phi) / pi;
      return 0.2 + 0.3 * t * t * t * (10 - 15 * t + 6 * t * t);
    };
    const int n = 2000;
    double removed = 0;
    for (int i = 0; i <= n; ++i) {
      const double w = (i == 0 || i == n) ? 1 : (i % 2 ? 4 : 2);
      removed += w * fillet_ring(2.0, rr(2 * pi * i / n), true) / (2 * pi);
    }
    removed *= 2 * pi / n / 3;
    CHECK(std::abs(f.volume() - (cyl.volume() - removed)) <= 1e-6);
    EdgeBlendOptions bad;
    bad.sizes = {{{2, 2, 2}, 0.3}};  // not on the edge
    CHECK_THROWS_AS(fillet_edges(box, e, r0, bad), GeometryError);
    bad.sizes = {{{0, 0, 2}, -1.0}};
    CHECK_THROWS_AS(chamfer_edges(box, e, r0, bad), GeometryError);
    // radii of a three-edge corner must agree
    const auto corner = edges_of(box, {"top", "front", "top", "left", "front", "left"});
    EdgeBlendOptions mix;
    mix.sizes = {
        {{0, 0, 2}, 0.3}, {{4, 0, 2}, 0.3}, {{0, 3, 2}, 0.3}, {{0, 0, 0}, 0.3}, {{0, 1.5, 2}, 0.4}};
    const Shell g = fillet_edges(box, corner, 0.3, mix);  // 0.3 at the corner, 0.4 on top-left
    CHECK(g.is_closed());
    mix.sizes.push_back({{0, 0, 1.9}, 0.5});
    CHECK_THROWS_AS(fillet_edges(box, corner, 0.3, mix), GeometryError);
  }
}
