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
