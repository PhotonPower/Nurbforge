#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <numbers>
#include <vector>

#include "nurbforge/core/error.hpp"
#include "nurbforge/geometry/primitives.hpp"
#include "nurbforge/modeling/construct.hpp"
#include "nurbforge/topology/shell.hpp"

using namespace nf;
constexpr double pi = std::numbers::pi;

namespace {

bool rel_close(double a, double b, double rel) { return std::abs(a - b) <= rel * std::abs(b); }

NurbsSurface sphere(double R) {
  return revolve(make_arc({0, 0, 0}, {0, 0, -1}, {1, 0, 0}, R, 0.0, pi), {0, 0, 0}, {0, 0, 1}, 2 * pi);
}

std::vector<TrimmedSurface> cube_faces(bool skip_top = false) {
  // unit cube, faces with mixed orientations on purpose
  std::vector<TrimmedSurface> f{
      TrimmedSurface(make_rectangle({0, 0, 0}, {1, 0, 0}, {0, 1, 0})),   // bottom
      TrimmedSurface(make_rectangle({0, 0, 0}, {1, 0, 0}, {0, 0, 1})),   // front
      TrimmedSurface(make_rectangle({0, 1, 0}, {0, 0, 1}, {1, 0, 0})),   // back
      TrimmedSurface(make_rectangle({0, 0, 0}, {0, 1, 0}, {0, 0, 1})),   // left
      TrimmedSurface(make_rectangle({1, 0, 0}, {0, 0, 1}, {0, 1, 0}))};  // right
  if (!skip_top) f.emplace_back(make_rectangle({0, 0, 1}, {0, 1, 0}, {1, 0, 0}));  // top
  return f;
}

}  // namespace

TEST_CASE("full sphere: seam edge and poles", "[shell]") {
  const double R = 1.5;
  const Shell s = Shell::sew({TrimmedSurface(sphere(R))}).oriented();
  CHECK(s.is_closed());
  CHECK(s.is_consistently_oriented());
  REQUIRE(s.edges().size() == 1);  // the seam, used twice by the same face
  CHECK(s.edges()[0].uses.size() == 2);
  CHECK(rel_close(s.area(), 4 * pi * R * R, 1e-7));
  CHECK(rel_close(s.volume(), 4.0 / 3.0 * pi * R * R * R, 1e-7));
}

TEST_CASE("cube from six rectangles", "[shell]") {
  const Shell raw = Shell::sew(cube_faces());
  CHECK(raw.edges().size() == 12);
  CHECK(raw.vertices().size() == 8);
  CHECK(raw.is_closed());
  CHECK(raw.free_edges().empty());
  const Shell s = raw.oriented();
  CHECK(s.is_consistently_oriented());
  CHECK(std::abs(s.volume() - 1.0) < 1e-10);
  CHECK(std::abs(s.area() - 6.0) < 1e-10);
}

TEST_CASE("open box has free edges", "[shell]") {
  const Shell s = Shell::sew(cube_faces(true));
  CHECK_FALSE(s.is_closed());
  CHECK(s.free_edges().size() == 4);
  CHECK_THROWS_AS(s.volume(), GeometryError);
}

TEST_CASE("spherical cap closed by a disk", "[shell][ssi]") {
  const double R = 2.0, h = 0.8;
  const auto ball = sphere(R);
  const auto plane = make_rectangle({-3, -3, h}, {6, 0, 0}, {0, 6, 0});
  auto caps = split(ball, plane);
  auto plates = split(plane, ball);
  REQUIRE(caps.size() == 2);
  REQUIRE(plates.size() == 2);
  // the upper (smaller) cap and the disk
  const auto& cap = caps[0].area(1e-6) < caps[1].area(1e-6) ? caps[0] : caps[1];
  const auto& disk = plates[0].area(1e-6) < plates[1].area(1e-6) ? plates[0] : plates[1];
  const Shell s = Shell::sew({cap, disk}).oriented();
  CHECK(s.is_closed());
  CHECK(s.is_consistently_oriented());
  const double hc = R - h;
  CHECK(rel_close(s.volume(1e-8), pi * hc * hc * (3 * R - hc) / 3.0, 1e-6));
}

TEST_CASE("cylinder with two trimmed lids", "[shell]") {
  const double R = 1.0, H = 2.0;
  const auto side = extrude(make_circle({0, 0, 0}, {0, 0, 1}, R), {0, 0, H});
  const auto square = [&](double z) { return make_rectangle({-R, -R, z}, {2 * R, 0, 0}, {0, 2 * R, 0}); };
  const TrimmedSurface bottom(square(0.0), circle_loop(0.5, 0.5, 0.5));
  const TrimmedSurface top(square(H), circle_loop(0.5, 0.5, 0.5));
  const Shell s = Shell::sew({TrimmedSurface(side), bottom, top}).oriented();
  CHECK(s.is_closed());
  CHECK(s.is_consistently_oriented());
  CHECK(rel_close(s.volume(1e-9), pi * R * R * H, 1e-7));
  CHECK(rel_close(s.area(1e-9), 2 * pi * R * H + 2 * pi * R * R, 1e-7));
}

TEST_CASE("three faces on one edge are non-manifold", "[shell]") {
  const Shell s = Shell::sew({TrimmedSurface(make_rectangle({0, 0, 0}, {1, 0, 0}, {0, 1, 0})),
                              TrimmedSurface(make_rectangle({0, 0, 0}, {1, 0, 0}, {0, -1, 0})),
                              TrimmedSurface(make_rectangle({0, 0, 0}, {1, 0, 0}, {0, 0, 1}))});
  CHECK(s.non_manifold_edges().size() == 1);
  CHECK_FALSE(s.is_closed());
}
