#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <numbers>
#include <thread>
#include <vector>

#include "nurbforge/core/error.hpp"
#include "nurbforge/geometry/primitives.hpp"
#include "nurbforge/modeling/construct.hpp"
#include "nurbforge/topology/boolean.hpp"
#include "nurbforge/topology/solid_ops.hpp"
#include "nurbforge/topology/solids.hpp"

using namespace nf;
constexpr double pi = std::numbers::pi;

namespace {

bool rel_close(double a, double b, double rel) { return std::abs(a - b) <= rel * std::abs(b); }

void check_solid(const Shell& s, double volume) {
  CHECK(s.is_closed());
  CHECK(s.is_consistently_oriented());
  CHECK(rel_close(s.volume(1e-10), volume, 1e-8));
}

}  // namespace

TEST_CASE("cap: planar openings become trimmed planar faces", "[solid_ops][cap]") {
  SECTION("open cylinder: two discs") {
    const auto side = extrude(make_circle({0, 0, 0}, {0, 0, 1}, 1.5), {0, 0, 2});
    const Shell tube = Shell::sew({TrimmedSurface(side)});
    REQUIRE(!tube.is_closed());
    const Shell solid = cap(tube);
    check_solid(solid, pi * 1.5 * 1.5 * 2);
    CHECK(solid.faces_labeled("cap").size() == 2);
    for (const auto i : solid.faces_labeled("cap"))
      CHECK(rel_close(solid.faces()[i].area(), pi * 2.25, 1e-9));
  }
  SECTION("pipe: annular caps with holes") {
    const auto outer = extrude(make_circle({0, 0, 0}, {0, 0, 1}, 2.0), {0, 0, 3});
    const auto inner = extrude(make_circle({0, 0, 0}, {0, 0, 1}, 1.0), {0, 0, 3});
    const Shell pipe = Shell::sew({TrimmedSurface(outer), TrimmedSurface(inner).flipped()});
    const Shell solid = cap(pipe, 1e-6, "end");
    check_solid(solid, pi * (4.0 - 1.0) * 3);
    REQUIRE(solid.faces_labeled("end").size() == 2);
    for (const auto i : solid.faces_labeled("end")) CHECK(solid.faces()[i].holes().size() == 1);
  }
  SECTION("box without its top and a hemisphere") {
    const Shell box = make_box({0, 0, 0}, {1, 2, 3});
    std::vector<TrimmedSurface> faces;
    for (const auto& f : box.faces())
      if (f.label() != "top") faces.push_back(f);
    check_solid(cap(Shell::sew(faces)), 6.0);
    const auto dome = revolve(make_arc({0, 0, 0}, {0, 0, 1}, {1, 0, 0}, 1.0, 0.0, pi / 2),
                              {0, 0, 0}, {0, 0, 1}, 2 * pi);
    check_solid(cap(Shell::sew({TrimmedSurface(dome)})), 2.0 / 3.0 * pi);
  }
  SECTION("closed shells are returned unchanged; non-planar openings are rejected") {
    const Shell ball = make_sphere({0, 0, 0}, 1.0);
    CHECK(cap(ball).faces().size() == ball.faces().size());
    const Shell saddle =
        Shell::sew({TrimmedSurface(make_bilinear({0, 0, 0}, {1, 0, 1}, {0, 1, 1}, {1, 1, 0}))});
    CHECK_THROWS_AS(cap(saddle), GeometryError);
  }
}

TEST_CASE("thicken: open shells and faces become solids with a wall thickness",
          "[solid_ops][thicken]") {
  SECTION("flat plate: a box, on either side and symmetric") {
    const Shell plate =
        Shell::sew({TrimmedSurface(make_rectangle({0, 0, 0}, {2, 0, 0}, {0, 3, 0}))});
    for (const double t : {0.1, -0.1}) {
      const Shell solid = thicken(plate, t);
      check_solid(solid, 0.6);
      CHECK(solid.faces_labeled("wall").size() == 4);
    }
    const Shell mid = thicken(plate, 0.2, true);
    check_solid(mid, 1.2);
    double zmin = 1e300, zmax = -1e300;
    for (const auto& p : mid.vertices()) zmin = std::min(zmin, p.z), zmax = std::max(zmax, p.z);
    CHECK(std::abs(zmin + 0.1) <= 1e-9);
    CHECK(std::abs(zmax - 0.1) <= 1e-9);
  }
  SECTION("cylinder wall: a tube; trimmed face: a disc with a hole") {
    const auto side = extrude(make_circle({0, 0, 0}, {0, 0, 1}, 1.0), {0, 0, 2});
    const Shell band = Shell::sew({TrimmedSurface(side)});
    const Shell tube = thicken(band, 0.1);  // along S_u x S_v, inward or outward
    const double out = pi * (1.21 - 1.0) * 2, in = pi * (1.0 - 0.81) * 2;
    CHECK(tube.is_closed());
    const double v = tube.volume(1e-10);
    CHECK((rel_close(v, out, 1e-6) || rel_close(v, in, 1e-6)));
    const auto sq = make_rectangle({-2, -2, 0}, {4, 0, 0}, {0, 4, 0});
    const Shell washer =
        Shell::sew({TrimmedSurface(sq, circle_loop(0.5, 0.5, 0.5), {circle_loop(0.5, 0.5, 0.25)})});
    check_solid(thicken(washer, 0.5), pi * (4.0 - 1.0) * 0.5);
  }
  SECTION("hemisphere: a spherical shell") {
    const auto dome = revolve(make_arc({0, 0, 0}, {0, 0, 1}, {1, 0, 0}, 1.0, 0.0, pi / 2),
                              {0, 0, 0}, {0, 0, 1}, 2 * pi);
    const Shell s = thicken(Shell::sew({TrimmedSurface(dome)}), 0.2);
    CHECK(s.is_closed());
    const double v = s.volume(1e-10);
    const double outward = 2.0 / 3.0 * pi * (1.2 * 1.2 * 1.2 - 1.0),
                 inward = 2.0 / 3.0 * pi * (1.0 - 0.8 * 0.8 * 0.8);
    CHECK((rel_close(v, outward, 1e-6) || rel_close(v, inward, 1e-6)));
  }
  SECTION("errors: closed shells, creases, zero thickness") {
    CHECK_THROWS_AS(thicken(make_sphere({0, 0, 0}, 1.0), 0.1), GeometryError);
    const Shell box = make_box({0, 0, 0}, {1, 1, 1});
    std::vector<TrimmedSurface> faces;
    for (const auto& f : box.faces())
      if (f.label() != "top") faces.push_back(f);
    CHECK_THROWS_AS(thicken(Shell::sew(faces), 0.1), GeometryError);  // creased open box
    const Shell plate =
        Shell::sew({TrimmedSurface(make_rectangle({0, 0, 0}, {1, 0, 0}, {0, 1, 0}))});
    CHECK_THROWS_AS(thicken(plate, 0.0), GeometryError);
  }
}

TEST_CASE("make_pipe: solid pipes along paths", "[solid_ops][pipe]") {
  SECTION("straight: cylinder and cone frustum") {
    const auto line = make_line({0, 0, 0}, {0, 0, 4});
    check_solid(make_pipe(line, 0.5), pi * 0.25 * 4);
    const Shell cone = make_pipe(line, 1.0, 0.5);
    check_solid(cone, pi * 4 / 3.0 * (1.0 + 0.5 + 0.25));
    REQUIRE(cone.faces_labeled("inlet").size() == 1);
    REQUIRE(cone.faces_labeled("outlet").size() == 1);
    CHECK(std::abs(cone.faces()[cone.faces_labeled("inlet")[0]].area() - pi) <= 1e-8);
    CHECK(std::abs(cone.faces()[cone.faces_labeled("outlet")[0]].area() - pi * 0.25) <= 1e-8);
  }
  SECTION("bent: quarter arc (Pappus) and a closed ring (torus)") {
    const double R = 3.0, r = 0.4;
    const auto bend = make_arc({0, 0, 0}, {1, 0, 0}, {0, 1, 0}, R, 0.0, pi / 2);
    const Shell elbow = make_pipe(bend, r);
    CHECK(elbow.is_closed());
    CHECK(rel_close(elbow.volume(1e-10), pi * r * r * R * pi / 2, 1e-6));
    const Shell ring = make_pipe(make_circle({0, 0, 0}, {0, 0, 1}, R), r);
    CHECK(ring.is_closed());
    CHECK(ring.faces_labeled("inlet").empty());
    CHECK(rel_close(ring.volume(1e-10), 2 * pi * pi * R * r * r, 1e-6));
  }
  SECTION("variable radius from a function") {
    const auto line = make_line({0, 0, 0}, {4, 0, 0});
    // r(t) = 0.5 + 0.2 sin(pi t): volume = pi int_0^4 r^2 dx
    const Shell bulge = make_pipe(line, [](double t) { return 0.5 + 0.2 * std::sin(pi * t); });
    // int_0^1 (0.5 + 0.2 sin(pi t))^2 dt = 0.25 + 2 * 0.5 * 0.2 * 2 / pi + 0.04 / 2
    const double I = 0.25 + 2 * 0.5 * 0.2 * (2 / pi) + 0.04 * 0.5;
    CHECK(rel_close(bulge.volume(1e-10), pi * 4 * I, 1e-6));
    CHECK_THROWS_AS(make_pipe(line, [](double t) { return 0.5 - t; }), GeometryError);
  }
}

TEST_CASE("remove_holes: small through and blind holes are filled", "[solid_ops][defeature]") {
  const Shell plate = make_box({-3, -2, 0}, {6, 4, 1});
  Shell part =
      boolean(plate, make_cylinder({-1.5, 0, -1}, {0, 0, 1}, 0.2, 3), BooleanOp::Difference);
  part = boolean(part, make_cylinder({1.5, 0, 0.5}, {0, 0, 1}, 0.25, 2), BooleanOp::Difference);
  part = boolean(part, make_cylinder({0, 0, -1}, {0, 0, 1}, 1.0, 3), BooleanOp::Difference);
  const double big = pi * 1.0 * 1.0, small = pi * 0.04, blind = pi * 0.0625 * 0.5;
  REQUIRE(rel_close(part.volume(1e-10), 24.0 - big - small - blind, 1e-6));
  // holes up to 0.6 across go, the bore of diameter 2 stays
  const Shell clean = remove_holes(part, 0.8);
  check_solid(clean, 24.0 - big);
  CHECK(clean.faces().size() < part.faces().size());
  // nothing to remove
  CHECK(remove_holes(part, 0.1).faces().size() == part.faces().size());
  CHECK_THROWS_AS(
      remove_holes(Shell::sew({TrimmedSurface(make_rectangle({0, 0, 0}, {1, 0, 0}, {0, 1, 0}))}),
                   1.0),
      GeometryError);
}

TEST_CASE("hollow: smooth solids become shells with a wall", "[solid_ops][hollow]") {
  // closed: a sphere with a concentric void
  const Shell ball = make_sphere({0, 0, 0}, 1.0);
  const Shell hollow_ball = hollow(ball, 0.1);
  check_solid(hollow_ball, 4.0 / 3.0 * pi * (1.0 - 0.729));
  CHECK(hollow_ball.components().size() == 2);
  // open: a sphere made of two labelled hemispheres, opened at the top -> a cup
  const auto upper = revolve(make_arc({0, 0, 0}, {0, 0, 1}, {1, 0, 0}, 1.0, 0.0, pi / 2), {0, 0, 0},
                             {0, 0, 1}, 2 * pi);
  const auto lower = revolve(make_arc({0, 0, 0}, {1, 0, 0}, {0, 0, -1}, 1.0, 0.0, pi / 2),
                             {0, 0, 0}, {0, 0, 1}, 2 * pi);
  const Shell two = Shell::sew({TrimmedSurface(upper).with_label("top"),
                                TrimmedSurface(lower).with_label("bottom")})
                        .oriented();
  REQUIRE(two.is_closed());
  const Shell cup = hollow(two, 0.1, {"top"});
  check_solid(cup, 2.0 / 3.0 * pi * (1.0 - 0.729));
  CHECK_THROWS_AS(hollow(two, 0.1, {"nothing"}), GeometryError);
  CHECK_THROWS_AS(hollow(make_box({0, 0, 0}, {1, 1, 1}), 0.1), GeometryError);  // creased
}

TEST_CASE("merge_faces: coplanar and same-surface pieces become one face", "[solid_ops][merge]") {
  SECTION("union of two overlapping boxes: 14 faces -> 10") {
    const Shell u =
        boolean(make_box({0, 0, 0}, {2, 2, 1}), make_box({1, 1, 0}, {2, 2, 1}), BooleanOp::Union);
    REQUIRE(u.faces().size() == 14);
    const Shell m = merge_faces(u);
    check_solid(m, 7.0);
    CHECK(m.faces().size() == 10);
    REQUIRE(m.faces_labeled("top").size() == 1);
    CHECK(rel_close(m.faces()[m.faces_labeled("top")[0]].area(), 7.0, 1e-12));
  }
  SECTION("two boxes side by side: the touching faces vanish, the rest merges into 6") {
    const Shell u =
        boolean(make_box({0, 0, 0}, {1, 1, 1}), make_box({1, 0, 0}, {1, 1, 1}), BooleanOp::Union);
    const Shell m = merge_faces(u);
    check_solid(m, 2.0);
    CHECK(m.faces().size() == 6);
  }
  SECTION("a plate with a bore keeps its hole; nothing to merge leaves the shell as it is") {
    const Shell part = boolean(make_box({-2, -1, 0}, {4, 2, 1}),
                               make_cylinder({0, 0, -1}, {0, 0, 1}, 0.4, 3), BooleanOp::Difference);
    const Shell m = merge_faces(part);
    CHECK(m.faces().size() == part.faces().size());
    check_solid(m, 8.0 - pi * 0.16);
  }
  SECTION("faces with different labels are not merged") {
    const Shell a = make_box({0, 0, 0}, {1, 1, 1});
    const Shell b =
        make_box({1, 0, 0}, {1, 1, 1})
            .with_face_labels({"bottom2", "top2", "front2", "back2", "left2", "right2"});
    const Shell m = merge_faces(boolean(a, b, BooleanOp::Union));
    check_solid(m, 2.0);
    CHECK(m.faces().size() == 10);
  }
}

TEST_CASE("constructive operations run concurrently on shared inputs", "[solid_ops][threads]") {
  // The Python bindings release the GIL for Booleans, splits, offsets, sweeps and the solid
  // operations, so they may run in parallel threads on the same (immutable) input objects.
  const Shell plate = make_box({-3, -2, 0}, {6, 4, 1});
  const Shell drill = make_cylinder({0, 0, -1}, {0, 0, 1}, 0.4, 3);
  const auto line = make_line({0, 0, 0}, {0, 0, 2});
  const auto run = [&] {
    const Shell part = boolean(plate, drill, BooleanOp::Difference);
    const auto halves = split(part, Plane{{0, 0, 0}, {1, 0, 0}});
    const Shell pipe = make_pipe(line, 0.3, 0.2);
    const Shell merged =
        merge_faces(boolean(plate, make_box({2, 1, 0}, {2, 2, 1}), BooleanOp::Union));
    return std::vector<double>{part.volume(1e-9), halves.first.volume(1e-9), pipe.volume(1e-9),
                               merged.volume(1e-9)};
  };
  const auto ref = run();
  std::vector<std::vector<double>> got(4);
  std::vector<std::thread> th;
  for (std::size_t k = 0; k < got.size(); ++k) th.emplace_back([&, k] { got[k] = run(); });
  for (auto& t : th) t.join();
  for (const auto& g : got) CHECK(g == ref);
}

TEST_CASE("slab along an open and a closed curve", "[solid_ops]") {
  // L-shaped wall, 0.5 wide and 1.5 high: strip area 1.25
  const Shell wall = make_slab(make_polyline({{0, 0, 0}, {2, 0, 0}, {2, 1, 0}}), 0.5, {0, 0, 1.5});
  CHECK(wall.is_closed());
  CHECK(std::abs(wall.volume() - 1.25 * 1.5) < 1e-9);
  CHECK(wall.faces_labeled("cap").size() == 2);
  // to the right of the curve: the mirror image at the corner
  const Shell right =
      make_slab(make_polyline({{0, 0, 0}, {2, 0, 0}, {2, 1, 0}}), -0.5, {0, 0, 1.5});
  CHECK(std::abs(right.volume() - 1.75 * 1.5) < 1e-9);
  // ring inside a circle of radius 1 (left of a counter-clockwise curve), wall 0.2
  const Shell ring = make_slab(make_circle({0, 0, 0}, {0, 0, 1}, 1.0), 0.2, {0, 0, 2});
  CHECK(ring.is_closed());
  CHECK(std::abs(ring.volume() - std::numbers::pi * (1.0 - 0.64) * 2.0) < 1e-5);
  CHECK_THROWS_AS(make_slab(make_circle({0, 0, 0}, {0, 0, 1}, 1.0), 0.0, {0, 0, 1}), GeometryError);
  CHECK_THROWS_AS(make_slab(make_circle({0, 0, 0}, {0, 0, 1}, 1.0), 0.2, {0, 0, 0}), GeometryError);
}

TEST_CASE("wire cut of a box", "[solid_ops]") {
  const Shell box = make_box({0, 0, 0}, {4, 2, 1});
  // a straight wire along y at x = 1.5, running in +y: the left part is x < 1.5
  const auto halves = wire_cut(box, make_line({1.5, -1, 0}, {1.5, 3, 0}), {0, 0, 1});
  REQUIRE(halves.size() == 2);
  CHECK(halves[0].is_closed());
  CHECK(std::abs(halves[0].volume() - 3.0) < 1e-9);
  CHECK(std::abs(halves[1].volume() - 5.0) < 1e-9);
  CHECK(halves[0].faces_labeled("cut").size() == 1);
  CHECK(halves[1].faces_labeled("cut").size() == 1);
  // a circular wire: the plug and the plate with the hole
  const auto plug = wire_cut(box, make_circle({2, 1, 5}, {0, 0, 1}, 0.5), {0, 0, 1});
  REQUIRE(plug.size() == 2);
  CHECK(std::abs(plug[0].volume() - std::numbers::pi * 0.25) < 1e-6);
  CHECK(std::abs(plug[1].volume() - (8.0 - std::numbers::pi * 0.25)) < 1e-6);
  // a curved wire: an arc of radius 5 about (2, -4), bulging to +y; both parts together are
  // the box
  const auto bent = wire_cut(
      box, make_arc({2, -4, 0}, {1, 0, 0}, {0, 1, 0}, 5.0, 0.5, std::numbers::pi - 0.5), {0, 0, 1});
  REQUIRE(bent.size() == 2);
  CHECK(std::abs(bent[0].volume() + bent[1].volume() - 8.0) < 1e-6);
  CHECK(bent[0].volume() > 0.5);
  CHECK(bent[1].volume() > 0.5);
  // a wire beside the solid cuts nothing
  const auto miss = wire_cut(box, make_line({-2, -1, 0}, {-2, 3, 0}), {0, 0, 1});
  REQUIRE(miss.size() == 1);
  CHECK(std::abs(miss[0].volume() - 8.0) < 1e-9);
}
