#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <numbers>

#include "nurbforge/analysis/drawing.hpp"
#include "nurbforge/analysis/properties.hpp"
#include "nurbforge/core/error.hpp"
#include "nurbforge/geometry/primitives.hpp"
#include "nurbforge/topology/boolean.hpp"
#include "nurbforge/topology/solids.hpp"

using namespace nf;
constexpr double pi = std::numbers::pi;

namespace {

double total_length(const std::vector<DrawingCurve>& cs) {
  double l = 0.0;
  for (const auto& c : cs) l += arc_length(c.curve);
  return l;
}

}  // namespace

TEST_CASE("visibility of a line behind a box", "[drawing]") {
  const Shell box = make_box({0, 0, 0}, {1, 1, 1});
  const NurbsCurve below = make_line({-1, 0.5, -1}, {2, 0.5, -1});
  const NurbsCurve above = make_line({-1, 0.5, 2}, {2, 0.5, 2});
  const auto pieces = visibility({below, above}, View::parallel({0, 0, -1}), {box});
  REQUIRE(pieces.size() == 4);
  CHECK(pieces[0].visible);
  CHECK(!pieces[1].visible);
  CHECK(pieces[2].visible);
  CHECK(std::abs(pieces[1].t0 - 1.0 / 3.0) < 1e-7);
  CHECK(std::abs(pieces[1].t1 - 2.0 / 3.0) < 1e-7);
  CHECK(pieces[1].source == 0);
  CHECK(distance(pieces[1].curve.evaluate(pieces[1].t0), {0, 0.5, -1}) < 1e-6);
  CHECK(pieces[3].visible);
  CHECK(pieces[3].source == 1);
  CHECK(pieces[3].t1 - pieces[3].t0 == 1.0);

  // from an eye point above the box the hidden part is wider: x in [-1/9, 1 + 1/9]
  const auto persp = visibility({below}, View::perspective({0.5, 0.5, 10}), {box});
  REQUIRE(persp.size() == 3);
  CHECK(!persp[1].visible);
  CHECK(std::abs(persp[1].t0 - (1.0 - 1.0 / 9.0) / 3.0) < 1e-6);
  CHECK(std::abs(persp[1].t1 - (2.0 + 1.0 / 9.0) / 3.0) < 1e-6);
}

TEST_CASE("edges of a box: three are hidden in a general view", "[drawing]") {
  const Shell box = make_box({0, 0, 0}, {1, 2, 3});
  std::vector<NurbsCurve> edges;
  for (const auto& e : box.edges()) edges.push_back(e.curve);
  REQUIRE(edges.size() == 12);
  const auto pieces = visibility(edges, View::parallel({-1, -2.2, -1.3}), {box});
  REQUIRE(pieces.size() == 12);  // no edge changes its visibility
  std::size_t hidden = 0;
  for (const auto& p : pieces)
    if (!p.visible) {
      ++hidden;
      // the hidden edges meet at the far corner, the origin
      CHECK(std::min(norm(p.curve.evaluate(p.t0)), norm(p.curve.evaluate(p.t1))) < 1e-9);
    }
  CHECK(hidden == 3);
}

TEST_CASE("isometric view of a cube: the hidden edges lie exactly behind visible ones",
          "[drawing]") {
  const Shell cube = make_box({0, 0, 0}, {1, 1, 1});
  std::vector<NurbsCurve> edges;
  for (const auto& e : cube.edges()) edges.push_back(e.curve);
  const auto pieces = visibility(edges, View::parallel({-1, -1, -1}), {cube});
  REQUIRE(pieces.size() == 12);
  CHECK(std::count_if(pieces.begin(), pieces.end(), [](const auto& p) { return !p.visible; }) == 3);
  const Drawing d = make2d({cube}, {-1, -1, -1});
  CHECK(d.visible.size() == 9);
  CHECK(d.hidden.size() == 3);  // the three spokes from the far corner
  const Outline hex = outline({cube}, {-1, -1, -1});
  REQUIRE(hex.loops.size() == 1);
  CHECK(hex.loops[0].size() == 6);
  CHECK(std::abs(region_area(hex.loops) - std::sqrt(3.0)) < 1e-9);
}

TEST_CASE("make2d of a box", "[drawing]") {
  const Shell box = make_box({0, 0, 0}, {1, 2, 3});
  const Drawing general = make2d({box}, {-1, -2.2, -1.3});
  CHECK(general.visible.size() == 9);
  CHECK(general.hidden.size() == 3);
  for (const auto& c : general.visible) {
    CHECK(c.kind == DrawingCurveKind::Edge);
    CHECK(std::abs(c.curve.evaluate(c.curve.domain_start()).z) < 1e-12);
    // the plane curve is the projection of the 3D piece
    const double t = 0.5 * (c.curve.domain_start() + c.curve.domain_end());
    const Vec3 q = c.curve.evaluate(t);
    const Vec3 p = general.origin + q.x * general.s_axis + q.y * general.t_axis;
    const Vec3 d = c.curve3d.evaluate(t) - p;
    CHECK(norm(cross(d, general.normal)) < 1e-9);
  }
  CHECK(std::abs(dot(general.s_axis, general.normal)) < 1e-12);
  CHECK(distance(cross(general.s_axis, general.t_axis), general.normal) < 1e-12);

  // top view: the outline rectangle only; the bottom edges lie under the top ones, the
  // vertical edges are points
  DrawingOptions opt;
  opt.s_axis = {1, 0, 0};
  const Drawing top = make2d({box}, {0, 0, -1}, opt);
  CHECK(top.visible.size() == 4);
  CHECK(top.hidden.empty());
  CHECK(std::abs(total_length(top.visible) - 6.0) < 1e-9);
  CHECK(distance(top.t_axis, {0, 1, 0}) < 1e-12);
  for (const auto& c : top.visible) {
    const Vec3 q = c.curve.evaluate(c.curve.domain_start());
    CHECK(q.x > -1e-9);
    CHECK(q.x < 1 + 1e-9);
    CHECK(q.y > -1e-9);
    CHECK(q.y < 2 + 1e-9);
  }
  CHECK_THROWS_AS(make2d({box}, {0, 0, 0}), GeometryError);
  opt.s_axis = {0, 0, 1};
  CHECK_THROWS_AS(make2d({box}, {0, 0, -1}, opt), GeometryError);
}

TEST_CASE("make2d of a cylinder from the side is a rectangle", "[drawing]") {
  const double R = 1.0, H = 3.0;
  const Shell cyl = make_cylinder({0, 0, 0}, {0, 0, 1}, R, H);
  const Drawing d = make2d({cyl}, {1, 0, 0});
  CHECK(d.hidden.empty());
  // two silhouette lines and the two rims seen edge-on (their far halves lie under the near)
  CHECK(std::abs(total_length(d.visible) - (2 * H + 4 * R)) < 1e-5);
  std::size_t silhouettes = 0;
  for (const auto& c : d.visible) {
    if (c.kind == DrawingCurveKind::Silhouette) ++silhouettes;
    for (int i = 0; i <= 8; ++i) {
      const Vec3 q = c.curve.evaluate(c.curve.domain_start() +
                                      (c.curve.domain_end() - c.curve.domain_start()) * i / 8.0);
      CHECK(std::abs(q.x) < R + 1e-6);
      CHECK(q.y > -1e-6);
      CHECK(q.y < H + 1e-6);
    }
  }
  CHECK(silhouettes == 2);
  // with the seam, a smooth edge
  DrawingOptions opt;
  opt.tangent_edges = true;
  const Drawing e = make2d({cyl}, {1, 0, 0}, opt);
  CHECK(std::abs(total_length(e.visible) + total_length(e.hidden) - (3 * H + 4 * R)) < 1e-5);
}

TEST_CASE("make2d: a box behind another one is hidden", "[drawing]") {
  const Shell front = make_box({0, 0, 2}, {2, 2, 1}), back = make_box({0.5, 0.5, 0}, {1, 1, 1});
  const Drawing d = make2d({front, back}, {0, 0, -1});
  CHECK(std::abs(total_length(d.visible) - 8.0) < 1e-9);
  CHECK(std::abs(total_length(d.hidden) - 4.0) < 1e-9);
  for (const auto& c : d.hidden) CHECK(c.shell == 1);
  DrawingOptions opt;
  opt.hidden = false;
  CHECK(make2d({front, back}, {0, 0, -1}, opt).hidden.empty());
  // partly hidden: the back box sticks out to the right
  const Shell wide = make_box({1, 0.5, 0}, {2, 1, 1});
  const Drawing e = make2d({front, wide}, {0, 0, -1});
  CHECK(std::abs(total_length(e.visible) - (8.0 + 3.0)) < 1e-6);
  CHECK(std::abs(total_length(e.hidden) - 3.0) < 1e-6);
}

TEST_CASE("outline of solids", "[drawing]") {
  const Shell box = make_box({0, 0, 0}, {1, 2, 3});
  const Vec3 dir = normalized(Vec3{-1, -2.2, -1.3});
  const Outline hex = outline({box}, dir);
  REQUIRE(hex.loops.size() == 1);
  CHECK(hex.loops[0].size() == 6);
  CHECK(std::abs(region_area(hex.loops) -
                 (std::abs(dir.x) * 6 + std::abs(dir.y) * 3 + std::abs(dir.z) * 2)) < 1e-9);

  const Outline disc = outline({make_sphere({1, 2, 3}, 2.0)}, {1, 1, 0.5});
  REQUIRE(disc.loops.size() == 1);
  CHECK(std::abs(region_area(disc.loops) - 4 * pi) < 1e-4);

  const Outline rect = outline({make_cylinder({0, 0, 0}, {0, 0, 1}, 1.0, 3.0)}, {1, 0, 0});
  REQUIRE(rect.loops.size() == 1);
  CHECK(std::abs(region_area(rect.loops) - 6.0) < 1e-5);

  // a torus along its axis: an annulus (outer loop and hole)
  const Outline ring = outline({make_torus({0, 0, 0}, {0, 0, 1}, 3.0, 1.0)}, {0, 0, -1});
  REQUIRE(ring.loops.size() == 2);
  CHECK(std::abs(region_area(ring.loops) - pi * (16 - 4)) < 1e-4);

  // two overlapping boxes: the union of their projections
  const Outline both =
      outline({make_box({0, 0, 0}, {2, 2, 1}), make_box({1, 1, 3}, {2, 2, 1})}, {0, 0, -1});
  REQUIRE(both.loops.size() == 1);
  CHECK(std::abs(region_area(both.loops) - 7.0) < 1e-9);
}

TEST_CASE("drawing and outline of a plate with a through hole", "[drawing]") {
  // 4 x 2 x 1 plate with a hole of radius 0.5 along z
  const Shell plate = boolean(make_box({0, 0, 0}, {4, 2, 1}),
                              make_cylinder({2, 1, -1}, {0, 0, 1}, 0.5, 3), BooleanOp::Difference);
  DrawingOptions opt;
  opt.s_axis = {1, 0, 0};
  const Drawing top = make2d({plate}, {0, 0, -1}, opt);
  CHECK(top.hidden.empty());
  CHECK(std::abs(total_length(top.visible) - (12.0 + pi)) < 1e-5);  // rectangle and circle

  const Outline contour = outline({plate}, {0, 0, -1}, opt);
  REQUIRE(contour.loops.size() == 2);
  CHECK(std::abs(region_area(contour.loops) - (8.0 - pi * 0.25)) < 1e-6);

  // front view: the hole shows as two hidden lines (its silhouettes) through the plate
  const Drawing front = make2d({plate}, {0, 1, 0}, opt);
  CHECK(std::abs(total_length(front.visible) - 10.0) < 1e-5);
  CHECK(std::abs(total_length(front.hidden) - 2.0) < 1e-5);
  for (const auto& c : front.hidden) {
    CHECK(c.kind == DrawingCurveKind::Silhouette);
    const Vec3 q = c.curve.evaluate(c.curve.domain_start());
    CHECK(std::abs(std::abs(q.x - 2.0) - 0.5) < 1e-6);
  }
  const Outline side = outline({plate}, {0, 1, 0}, opt);
  REQUIRE(side.loops.size() == 1);
  CHECK(std::abs(region_area(side.loops) - 4.0) < 1e-6);
}

TEST_CASE("perspective drawing and outline", "[drawing]") {
  // cube of edge 2 about the origin, eye on the z axis at distance 10, plane z = 0
  const Shell cube = make_box({-1, -1, -1}, {2, 2, 2});
  const View view = View::perspective({0, 0, 10});
  DrawingOptions opt;
  opt.s_axis = {1, 0, 0};
  const Drawing d = make2d({cube}, view, opt);
  const double near = 10.0 / 9.0, far = 10.0 / 11.0;  // half sizes of the projected faces
  CHECK(std::abs(total_length(d.visible) - 8 * near) < 1e-6);
  CHECK(std::abs(total_length(d.hidden) - (8 * far + 4 * std::sqrt(2.0) * (near - far))) < 1e-6);
  CHECK(distance(d.normal, {0, 0, 1}) < 1e-12);
  for (const auto& c : d.visible) {
    const Vec3 q = c.curve.evaluate(0.5 * (c.curve.domain_start() + c.curve.domain_end()));
    CHECK(std::abs(std::max(std::abs(q.x), std::abs(q.y)) - near) < 1e-9);
    // the plane curve is the central projection of the 3D piece
    const double t = c.curve.domain_start() + 0.3 * (c.curve.domain_end() - c.curve.domain_start());
    const Vec3 p = c.curve3d.evaluate(t), r = c.curve.evaluate(t);
    CHECK(std::abs(r.x - 10 * p.x / (10 - p.z)) < 1e-9);
    CHECK(std::abs(r.y - 10 * p.y / (10 - p.z)) < 1e-9);
  }
  const Outline square = outline({cube}, view, opt);
  REQUIRE(square.loops.size() == 1);
  CHECK(std::abs(region_area(square.loops) - 4 * near * near) < 1e-9);

  // a unit sphere from distance 5: a circle of radius 1 / sqrt(1 - 1/25)
  const Outline disc = outline({make_sphere({0, 0, 0}, 1.0)}, View::perspective({3, 0, 4}));
  REQUIRE(disc.loops.size() == 1);
  CHECK(std::abs(region_area(disc.loops) - pi / 0.96) < 1e-4);

  CHECK_THROWS_AS(make2d({cube}, View::perspective({0, 0, 0.5})), GeometryError);  // eye inside
  CHECK_THROWS_AS(make2d({cube}, View::perspective({0, 0, 0})), GeometryError);    // at the origin
}

TEST_CASE("visibility changes are found at the crossings with the occluder's outline",
          "[drawing]") {
  // a strip 1 mm wide above a line 3 m long: far narrower than the sample spacing
  const Shell strip = make_box({1.2345, -1, 1}, {0.001, 2, 0.1});
  const NurbsCurve line = make_line({0, 0, 0}, {3, 0, 0});
  VisibilityOptions opt;
  opt.samples = 8;
  const auto pieces = visibility({line}, View::parallel({0, 0, -1}), {strip}, opt);
  REQUIRE(pieces.size() == 3);
  CHECK(!pieces[1].visible);
  CHECK(std::abs(pieces[1].t0 - 1.2345 / 3) < 1e-9);
  CHECK(std::abs(pieces[1].t1 - 1.2355 / 3) < 1e-9);
  // sampling alone misses it
  opt.events = false;
  CHECK(visibility({line}, View::parallel({0, 0, -1}), {strip}, opt).size() == 1);
  // a line passing through a face changes its visibility away from any outline
  const Shell plate = make_box({-1, -1, 0}, {2, 2, 0.1});
  const NurbsCurve needle = make_line({0.2, 0.1, -1}, {0.3, 0.2, 1});
  const auto through = visibility({needle}, View::parallel({0, 0, -1}), {plate});
  REQUIRE(through.size() == 2);
  CHECK(!through[0].visible);
  CHECK(through[1].visible);
  CHECK(std::abs(needle.evaluate(through[0].t1).z - 0.1) < 1e-5);  // leaves through the top
}
