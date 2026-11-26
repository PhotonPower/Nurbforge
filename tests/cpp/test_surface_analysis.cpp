#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <numbers>

#include "nurbforge/analysis/surface_analysis.hpp"
#include "nurbforge/core/error.hpp"
#include "nurbforge/geometry/primitives.hpp"
#include "nurbforge/modeling/construct.hpp"
#include "nurbforge/modeling/surface_tools.hpp"
#include "nurbforge/topology/edge_blend.hpp"
#include "nurbforge/topology/solids.hpp"

using namespace nf;
constexpr double pi = std::numbers::pi;

TEST_CASE("curvature maps of sphere, cylinder and torus", "[surface_analysis]") {
  const NurbsSurface sphere = make_sphere({0, 0, 0}, 2.0).faces()[0].surface();
  const CurvatureMap s = curvature_map(sphere, 4);
  REQUIRE(s.gaussian.size() == s.u.size() * s.v.size());
  CHECK(std::abs(s.gaussian_range.min - 0.25) < 1e-6);
  CHECK(std::abs(s.gaussian_range.max - 0.25) < 1e-6);
  CHECK(std::abs(std::abs(s.mean_range.min) - 0.5) < 1e-6);
  CHECK(std::abs(s.min_radius - 2.0) < 1e-5);

  const NurbsSurface cyl = extrude(make_circle({0, 0, 0}, {0, 0, 1}, 1.5), {0, 0, 4});
  const CurvatureMap c = curvature_map(cyl);
  CHECK(std::abs(c.gaussian_range.min) < 1e-12);
  CHECK(std::abs(c.gaussian_range.max) < 1e-12);
  CHECK(std::abs(std::abs(c.mean_range.max) - 1.0 / 3.0) < 1e-12);
  CHECK(std::abs(c.min_radius - 1.5) < 1e-12);

  // torus R = 3, r = 1: K from -1 / (r (R - r)) on the inner to 1 / (r (R + r)) on the outer
  // equator
  const NurbsSurface torus = make_torus({0, 0, 0}, {0, 0, 1}, 3.0, 1.0).faces()[0].surface();
  const CurvatureMap t = curvature_map(torus, 4);
  CHECK(std::abs(t.gaussian_range.min + 0.5) < 1e-9);
  CHECK(std::abs(t.gaussian_range.max - 0.25) < 1e-9);
  const Vec3 inner = torus.evaluate(t.gaussian_range.min_u, t.gaussian_range.min_v);
  CHECK(std::abs(std::hypot(inner.x, inner.y) - 2.0) < 1e-9);
  CHECK(std::abs(t.min_radius - 1.0) < 1e-9);

  // trimmed: only the samples inside count
  const auto du = torus.domain_u(), dv = torus.domain_v();
  const Vec3 outer = torus.evaluate(t.gaussian_range.max_u, t.gaussian_range.max_v);
  CHECK(std::abs(std::hypot(outer.x, outer.y) - 4.0) < 1e-9);
  const double vo = t.gaussian_range.max_v, span = 0.1 * (dv[1] - dv[0]);
  const double v0 = std::clamp(vo - span, dv[0], dv[1] - 2 * span);
  const TrimmedSurface band(torus, rectangle_loop(du[0], du[1], v0, v0 + 2 * span));
  const CurvatureMap b = curvature_map(band, 8);
  CHECK(std::count(b.inside.begin(), b.inside.end(), true) > 0);
  CHECK(std::count(b.inside.begin(), b.inside.end(), false) > 0);
  CHECK(b.gaussian_range.min > 0.0);  // the outer band only
  CHECK_THROWS_AS(curvature_map(torus, 0), GeometryError);
}

TEST_CASE("draft angles", "[surface_analysis]") {
  // cone with its apex up: the wall makes the angle atan(R / h) with the axis
  const double R = 1.0, h = 2.0;
  const NurbsSurface cone = extrude_to_point(make_circle({0, 0, 0}, {0, 0, 1}, R), {0, 0, h});
  const double expect = std::atan(R / h);
  const DraftAngleMap m = draft_angle_map(cone, {0, 0, 1}, 4);
  CHECK(std::abs(std::abs(m.range.min) - expect) < 1e-9);
  CHECK(std::abs(std::abs(m.range.max) - expect) < 1e-9);
  CHECK(m.range.min * m.range.max > 0.0);  // one sign everywhere
  CHECK(std::abs(draft_angle(cone, 0.3, 0.4, {0, 0, -5}) + draft_angle(cone, 0.3, 0.4, {0, 0, 1})) <
        1e-12);
  // a wall parallel to the pull direction has no draft; a sphere has every angle
  const NurbsSurface cyl = extrude(make_circle({0, 0, 0}, {0, 0, 1}, 1.0), {0, 0, 2});
  const DraftAngleMap w = draft_angle_map(cyl, {0, 0, 1});
  CHECK(std::abs(w.range.min) < 1e-12);
  CHECK(std::abs(w.range.max) < 1e-12);
  const DraftAngleMap s =
      draft_angle_map(make_sphere({0, 0, 0}, 1.0).faces()[0].surface(), {0, 0, 1}, 4);
  CHECK(std::abs(s.range.min + pi / 2) < 1e-6);
  CHECK(std::abs(s.range.max - pi / 2) < 1e-6);
  CHECK_THROWS_AS(draft_angle(cyl, 0.5, 0.5, {0, 0, 0}), GeometryError);
}

TEST_CASE("surface direction, closed sides and poles", "[surface_analysis]") {
  const NurbsSurface plane = make_rectangle({0, 0, 0}, {2, 0, 0}, {0, 3, 0});
  const SurfaceDirection p = surface_direction(plane);
  CHECK(distance(p.u_dir, {1, 0, 0}) < 1e-12);
  CHECK(distance(p.v_dir, {0, 1, 0}) < 1e-12);
  CHECK(distance(p.normal, {0, 0, 1}) < 1e-12);
  CHECK(distance(p.point, {1, 1.5, 0}) < 1e-12);
  CHECK(p.planar);
  CHECK(!p.closed_u);
  CHECK(!p.closed_v);
  CHECK(distance(surface_direction(plane.transposed()).normal, {0, 0, -1}) < 1e-12);

  const SurfaceDirection c =
      surface_direction(extrude(make_circle({0, 0, 0}, {0, 0, 1}, 1.0), {0, 0, 2}));
  CHECK(c.closed_u);
  CHECK(!c.closed_v);
  CHECK(!c.planar);
  CHECK(!(c.pole[0] || c.pole[1] || c.pole[2] || c.pole[3]));

  const NurbsSurface sphere = make_sphere({0, 0, 0}, 1.0).faces()[0].surface();
  const SurfaceDirection s = surface_direction(sphere);
  CHECK((s.closed_u != s.closed_v));  // closed around the axis only
  CHECK(std::count(std::begin(s.pole), std::end(s.pole), true) == 2);
}

TEST_CASE("continuity between surfaces", "[surface_analysis]") {
  // floor z = 0 for x <= 0, continued by a quarter cylinder of radius 1 about (0, y, 1)
  const NurbsSurface floor = make_rectangle({-2, 0, 0}, {2, 0, 0}, {0, 3, 0});
  const NurbsSurface bend =
      extrude(make_arc({0, 0, 1}, {0, 0, -1}, {1, 0, 0}, 1.0, 0.0, pi / 2), {0, 3, 0});
  const SurfaceContinuity g1 = continuity(floor, SurfaceSide::UMax, bend);
  CHECK(g1.geometric == 1);
  CHECK(g1.max_gap < 1e-12);
  CHECK(g1.max_angle < 1e-12);
  CHECK(std::abs(g1.max_curvature_difference - 1.0) < 1e-9);  // curvature 0 against 1

  const NurbsSurface wall = make_rectangle({0, 0, 0}, {0, 0, 2}, {0, 3, 0});
  const SurfaceContinuity g0 = continuity(floor, SurfaceSide::UMax, wall);
  CHECK(g0.geometric == 0);
  CHECK(std::abs(g0.max_angle - pi / 2) < 1e-12);
  CHECK(continuity(floor, SurfaceSide::UMin, wall).geometric == -1);

  // the two parts of one surface are curvature continuous
  const NurbsSurface patch = make_sphere({0, 0, 0}, 2.0).faces()[0].surface();
  const auto dv = patch.domain_v();
  const auto [low, high] = patch.split_v(dv[0] + 0.37 * (dv[1] - dv[0]));
  const SurfaceContinuity g2 = continuity(low, SurfaceSide::VMax, high);
  CHECK(g2.geometric == 2);
  CHECK(g2.max_curvature_difference < 1e-6);
  // also for surfaces of opposite orientation
  CHECK(continuity(low, SurfaceSide::VMax, high.transposed()).geometric == 2);
}

TEST_CASE("continuity across the edges of shells", "[surface_analysis]") {
  const Shell box = make_box({0, 0, 0}, {1, 2, 3});
  const auto b = edge_continuity(box);
  REQUIRE(b.size() == 12);
  for (const auto& k : b) {
    CHECK(k.geometric == 0);
    CHECK(std::abs(k.max_angle - pi / 2) < 1e-9);
  }
  // cylinder: the seam is smooth, the rims are sharp
  const Shell cyl = make_cylinder({0, 0, 0}, {0, 0, 1}, 1.0, 2.0);
  const auto c = edge_continuity(cyl);
  REQUIRE(c.size() == 3);
  CHECK(std::count_if(c.begin(), c.end(), [](const auto& k) { return k.geometric == 2; }) == 1);
  CHECK(std::count_if(c.begin(), c.end(), [](const auto& k) { return k.geometric == 0; }) == 2);
  // a filleted edge: the borders of the fillet are tangent edges
  const Shell round = fillet_edges(box, edges_between(box, "top", "front"), 0.3);
  const auto r = edge_continuity(round);
  CHECK(std::count_if(r.begin(), r.end(), [](const auto& k) { return k.geometric == 1; }) == 2);
  // an open shell: free edges have no continuity
  const Shell open = Shell::sew({TrimmedSurface(make_rectangle({0, 0, 0}, {1, 0, 0}, {0, 1, 0}))});
  for (const auto& k : edge_continuity(open)) CHECK(k.geometric == -1);
}

TEST_CASE("duplicate surfaces", "[surface_analysis]") {
  const NurbsSurface plane = make_rectangle({0, 0, 0}, {2, 0, 0}, {0, 3, 0});
  const NurbsSurface cyl = extrude(make_circle({0, 0, 0}, {0, 0, 1}, 1.0), {0, 0, 2});
  const std::vector<NurbsSurface> all{
      plane,                                               // 0
      cyl,                                                 // 1
      plane.transposed(),                                  // 2: the same rectangle, flipped
      cyl.elevate_degree_u().elevate_degree_v(),           // 3: other representation
      make_rectangle({0, 0, 0}, {1, 0, 0}, {0, 3, 0}),     // 4: half of the plane
      make_rectangle({0, 0, 1e-3}, {2, 0, 0}, {0, 3, 0}),  // 5: shifted
  };
  const auto groups = find_duplicate_surfaces(all, 1e-9);
  REQUIRE(groups.size() == 2);
  CHECK(groups[0] == std::vector<std::size_t>{0, 2});
  CHECK(groups[1] == std::vector<std::size_t>{1, 3});
  CHECK(find_duplicate_surfaces(all, 2e-3)[0] == std::vector<std::size_t>{0, 2, 5});
  CHECK_THROWS_AS(find_duplicate_surfaces(all, -1.0), GeometryError);
}
