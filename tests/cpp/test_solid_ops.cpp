#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <numbers>
#include <vector>

#include "nurbforge/core/error.hpp"
#include "nurbforge/geometry/primitives.hpp"
#include "nurbforge/modeling/construct.hpp"
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
