#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <numbers>

#include "nurbforge/analysis/closest_point.hpp"
#include "nurbforge/analysis/properties.hpp"
#include "nurbforge/core/error.hpp"
#include "nurbforge/fitting/fit.hpp"
#include "nurbforge/geometry/primitives.hpp"
#include "nurbforge/modeling/blend.hpp"
#include "nurbforge/modeling/construct.hpp"
#include "nurbforge/modeling/fillet.hpp"
#include "nurbforge/modeling/project.hpp"
#include "nurbforge/nurbs/join.hpp"

using namespace nf;
constexpr double pi = std::numbers::pi;

TEST_CASE("chamfers between curves and at corners", "[chamfer]") {
  SECTION("polyline corner: asymmetric distances, trimmed chain") {
    const auto a = make_line({0, 0, 0}, {2, 0, 0}), b = make_line({2, 0, 0}, {2, 2, 0});
    const auto c = chamfer_curves(a, b, 0.5, 0.25);
    CHECK(distance(c.point_a, Vec3{1.5, 0, 0}) <= 1e-14);
    CHECK(distance(c.point_b, Vec3{2, 0.25, 0}) <= 1e-14);
    CHECK(distance(c.line.evaluate(0.0), c.point_a) <= 1e-15);
    CHECK(distance(c.line.evaluate(1.0), c.point_b) <= 1e-15);
    REQUIRE(c.trimmed_a);
    REQUIRE(c.trimmed_b);
    const auto chain = join(join(*c.trimmed_a, c.line), *c.trimmed_b);
    CHECK(std::abs(arc_length(chain) - (1.5 + std::hypot(0.5, 0.25) + 1.75)) <= 1e-12);
  }
  SECTION("crossing lines: the kept parts are the longer sides") {
    const auto a = make_line({0, 0, 0}, {4, 0, 0}), b = make_line({1, -1, 0}, {1, 3, 0});
    const auto c = chamfer_curves(a, b, 1.0, 1.0);
    CHECK(distance(c.point_a, Vec3{2, 0, 0}) <= 1e-12);
    CHECK(distance(c.point_b, Vec3{1, 1, 0}) <= 1e-12);
  }
  SECTION("arc and line: chamfer points at the exact distances; 3D curves") {
    const auto arc = make_arc({0, 0, 0}, {1, 0, 0}, {0, 0, 1}, 2.0, 0.0, pi / 2);  // xz plane
    const auto line = make_line({0, 0, 2}, {0, 3, 2});
    const auto c = chamfer_curves(arc, line, 0.7, 0.4);
    const Vec3 X{0, 0, 2};
    CHECK(std::abs(distance(c.point_a, X) - 0.7) <= 1e-13);
    CHECK(std::abs(distance(c.point_b, X) - 0.4) <= 1e-13);
    CHECK(std::abs(norm(c.point_a) - 2.0) <= 1e-13);  // on the arc
  }
  SECTION("errors") {
    const auto a = make_line({0, 0, 0}, {1, 0, 0}), b = make_line({1, 0, 0}, {1, 1, 0});
    CHECK_THROWS_AS(chamfer_curves(a, b, 2.0, 0.5), GeometryError);
    CHECK_THROWS_AS(chamfer_curves(a, b, 0.0, 0.5), GeometryError);
    CHECK_THROWS_AS(chamfer_curves(a, make_line({0, 1, 0}, {1, 1, 0}), 0.1, 0.1), GeometryError);
  }
  SECTION("chamfer_corners on a closed rectangle") {
    const auto rect = make_rectangle_curve({0, 0, 0}, {4, 0, 0}, {0, 2, 0});
    const auto c = chamfer_corners(rect, 0.5);
    CHECK(distance(c.evaluate(c.domain_start()), c.evaluate(c.domain_end())) <= 1e-12);
    CHECK(std::abs(arc_length(c) - (12.0 - 8 * 0.5 + 4 * 0.5 * std::sqrt(2.0))) <= 1e-12);
    for (int i = 0; i <= 200; ++i) {  // never reaches a corner
      const Vec3 p = c.evaluate(c.domain_start() + (c.domain_end() - c.domain_start()) * i / 200);
      for (const Vec3 corner : {Vec3{0, 0, 0}, Vec3{4, 0, 0}, Vec3{4, 2, 0}, Vec3{0, 2, 0}})
        CHECK(distance(p, corner) >= 0.5 / std::sqrt(2.0) - 1e-12);
    }
    CHECK_THROWS_AS(chamfer_corners(rect, 1.5), GeometryError);  // overlaps on the short sides
  }
}

TEST_CASE("best-fit line and plane", "[fit]") {
  SECTION("line through noisy points") {
    const Vec3 p0{1, 2, 3}, d = normalized(Vec3{1, 1, 0}), e = Vec3{0, 0, 1};
    std::vector<Vec3> pts;
    for (int i = -5; i <= 5; ++i) pts.push_back(p0 + i * d + ((i % 2 == 0) ? 0.01 : -0.01) * e);
    const auto f = fit_line(pts);
    CHECK(norm(cross(f.direction, d)) <= 1e-3);
    CHECK(f.direction.x > 0.0);
    CHECK(std::abs(f.max_deviation - 0.01) <= 2e-3);
    CHECK_THROWS_AS(fit_line({p0, p0}), GeometryError);
  }
  SECTION("plane through exact points") {
    std::vector<Vec3> pts;
    for (int i = 0; i < 5; ++i)
      for (int j = 0; j < 4; ++j) {
        const double x = i * 0.7 - 1, y = j * 1.3 + 0.2;
        pts.push_back({x, y, 0.5 * x + 0.2 * y + 1.0});
      }
    const auto f = fit_plane(pts);
    CHECK(norm(cross(f.plane.normal, normalized(Vec3{-0.5, -0.2, 1}))) <= 1e-12);
    CHECK(f.plane.normal.z > 0.0);
    CHECK(f.max_deviation <= 1e-12);
    CHECK(f.rms <= 1e-12);
    CHECK_THROWS_AS(fit_plane({{0, 0, 0}, {1, 1, 1}, {2, 2, 2}}), GeometryError);
  }
}

TEST_CASE("parallel projection onto a plane", "[project]") {
  const auto circle = make_circle({1, 2, 3}, {0, 1, 1}, 2.0);
  const Plane xy{{0, 0, -1}, {0, 0, 1}};
  const auto ortho = project_to_plane(circle, xy);
  const Vec3 dir{1, 0, -2};
  const auto oblique = project_to_plane(circle, xy, dir);
  for (int i = 0; i <= 40; ++i) {
    const double t = i / 40.0;
    const Vec3 p = circle.evaluate(t), q = ortho.evaluate(t), r = oblique.evaluate(t);
    CHECK(distance(q, Vec3{p.x, p.y, -1}) <= 1e-13);
    CHECK(std::abs(r.z + 1) <= 1e-13);
    CHECK(norm(cross(r - p, dir)) <= 1e-12);
  }
  const auto patch = make_rectangle({0, 0, 1}, {1, 0, 1}, {0, 1, 0});
  const auto flat = project_to_plane(patch, xy);
  CHECK(std::abs(flat.evaluate(0.5, 0.5).z + 1) <= 1e-14);
  CHECK_THROWS_AS(project_to_plane(circle, xy, {1, 0, 0}), GeometryError);
}

TEST_CASE("tween curves and tangent constructions", "[blend]") {
  SECTION("tween between concentric circles is a circle") {
    const auto a = make_circle({0, 0, 0}, {0, 0, 1}, 1.0),
               b = make_circle({0, 0, 0}, {0, 0, 1}, 3.0);
    const auto tw = tween_curves(a, b, 3);
    REQUIRE(tw.size() == 3);
    for (int k = 0; k < 3; ++k)
      for (int i = 0; i <= 20; ++i)
        CHECK(std::abs(norm(tw[static_cast<std::size_t>(k)].evaluate(i / 20.0)) -
                       (1.5 + 0.5 * k)) <= 1e-12);
    CHECK_THROWS_AS(tween_curves(a, b, 0), GeometryError);
  }
  SECTION("tween between a line and a cubic") {
    const auto a = make_line({0, 0, 0}, {1, 0, 0});
    const auto b = NurbsCurve::from_points(3, {{0, 1, 0}, {0.3, 2, 0}, {0.6, 0, 0}, {1, 1, 0}});
    const auto tw = tween_curves(a, b, 1);
    for (int i = 0; i <= 10; ++i) {
      const double t = i / 10.0;
      CHECK(distance(tw[0].evaluate(t), 0.5 * (a.evaluate(t) + b.evaluate(t))) <= 1e-14);
    }
  }
  SECTION("arc from start, tangent and end") {
    const auto q = make_arc_tangent({0, 0, 0}, {1, 0, 0}, {1, 1, 0});
    CHECK(distance(q.evaluate(1.0), Vec3{1, 1, 0}) <= 1e-15);
    CHECK(distance(normalized(q.derivatives(0.0, 1)[1]), Vec3{1, 0, 0}) <= 1e-14);
    const auto big = make_arc_tangent({0, 0, 0}, {1, 0, 0}, {-1, 1, 0});  // 270 degrees
    for (int i = 0; i <= 20; ++i)
      CHECK(std::abs(distance(big.evaluate(i / 20.0), Vec3{0, 1, 0}) - 1.0) <= 1e-13);
    CHECK(std::abs(arc_length(big) - 1.5 * pi) <= 1e-9);
    CHECK_THROWS_AS(make_arc_tangent({0, 0, 0}, {1, 0, 0}, {2, 0, 0}), GeometryError);
  }
  SECTION("tangents from a point to a circle") {
    const auto circle = make_circle({0, 0, 0}, {0, 0, 1}, 1.0);
    const auto t = tangent_lines(Vec3{2, 0, 0}, circle);
    REQUIRE(t.size() == 2);
    for (const auto& l : t) {
      CHECK(std::abs(l.point_b.x - 0.5) <= 1e-12);
      CHECK(std::abs(std::abs(l.point_b.y) - std::sqrt(3.0) / 2) <= 1e-12);
      CHECK(std::isnan(l.param_a));
    }
    CHECK(tangent_lines(Vec3{0.2, 0, 0}, circle).empty());  // inside: no tangents
    CHECK_THROWS_AS(tangent_lines(Vec3{2, 0, 1}, circle), GeometryError);
  }
  SECTION("common tangents of two circles") {
    const auto c1 = make_circle({0, 0, 0}, {0, 0, 1}, 1.0),
               c2 = make_circle({4, 0, 0}, {0, 0, 1}, 1.0);
    const auto t = tangent_lines(c1, c2);
    REQUIRE(t.size() == 4);
    int outer = 0;
    for (const auto& l : t) {
      CHECK(std::abs(norm(l.point_a) - 1.0) <= 1e-12);
      CHECK(std::abs(distance(l.point_b, Vec3{4, 0, 0}) - 1.0) <= 1e-12);
      const Vec3 d = normalized(l.point_b - l.point_a);
      CHECK(std::abs(dot(d, l.point_a)) <= 1e-10);                  // perpendicular to radius 1
      CHECK(std::abs(dot(d, l.point_b - Vec3{4, 0, 0})) <= 1e-10);  // and radius 2
      if (std::abs(l.point_a.y - l.point_b.y) <= 1e-10) ++outer;
    }
    CHECK(outer == 2);
  }
}

TEST_CASE("curves on surfaces", "[project]") {
  const auto cyl = extrude(make_circle({0, 0, 0}, {0, 0, 1}, 1.0), {0, 0, 4});
  SECTION("interpolate through points on a cylinder") {
    std::vector<Vec3> pts;
    for (int i = 0; i <= 4; ++i) {
      const double a = 0.3 + 0.25 * i;
      pts.push_back({std::cos(a), std::sin(a), 0.5 + 0.6 * i * i / 4.0});
    }
    const auto r = interpolate_on_surface(cyl, pts, 3, false, 1e-8);
    for (const auto& p : pts) CHECK(closest_point(r.curve, p).distance <= 1e-8);
    CHECK(distance(r.curve.evaluate(0.0), pts.front()) <= 1e-8);
    CHECK(distance(r.curve.evaluate(1.0), pts.back()) <= 1e-8);
    for (int i = 0; i <= 50; ++i) {
      const Vec3 p = r.curve.evaluate(i / 50.0);
      CHECK(std::abs(std::hypot(p.x, p.y) - 1.0) <= 2e-8);
      const Vec3 w = r.pcurve.evaluate(i / 50.0);
      CHECK(distance(cyl.evaluate(w.x, w.y), p) <= 2e-8);
    }
  }
  SECTION("geodesic offset on a plane and a cylinder") {
    const auto plane = make_rectangle({0, 0, 0}, {4, 0, 0}, {0, 4, 0});
    const auto pc = make_line({0.1, 0.2, 0}, {0.9, 0.5, 0});
    const auto off = offset_on_surface(plane, pc, 0.3, 1e-9);
    const Vec3 a = plane.evaluate(0.1, 0.2), b = plane.evaluate(0.9, 0.5);
    const Vec3 left = normalized(cross(Vec3{0, 0, 1}, b - a));
    for (int i = 0; i <= 20; ++i) {
      const double t = i / 20.0;
      CHECK(distance(off.curve.evaluate(t), a + t * (b - a) + 0.3 * left) <= 1e-8);
    }
    // circle at z = 1 on the cylinder (u around): offset 0.5 moves it up to z = 1.5
    const auto ring = make_line({0, 0.25, 0}, {1, 0.25, 0});
    const auto up = offset_on_surface(cyl, ring, 0.5, 1e-8);
    const auto down = offset_on_surface(cyl, ring, -0.5, 1e-8);
    for (int i = 0; i <= 20; ++i) {
      const Vec3 p = up.curve.evaluate(i / 20.0), q = down.curve.evaluate(i / 20.0);
      CHECK(std::abs(p.z - 1.5) <= 1e-7);
      CHECK(std::abs(q.z - 0.5) <= 1e-7);
      CHECK(std::abs(std::hypot(p.x, p.y) - 1.0) <= 1e-7);
    }
    CHECK_THROWS_AS(offset_on_surface(cyl, ring, 5.0), GeometryError);  // leaves the surface
  }
  SECTION("geodesic offset of the equator on a sphere: a latitude circle") {
    const double R = 2.0, d = 0.8;
    const auto sphere = revolve(make_arc({0, 0, 0}, {1, 0, 0}, {0, 0, 1}, R, -pi / 2, pi / 2),
                                {0, 0, 0}, {0, 0, 1}, 2 * pi);
    const auto equator = make_line({0, 0.5, 0}, {1, 0.5, 0});
    const auto off = offset_on_surface(sphere, equator, d, 1e-8);
    for (int i = 0; i <= 20; ++i) {
      const Vec3 p = off.curve.evaluate(i / 20.0);
      CHECK(std::abs(p.z - R * std::sin(d / R)) <= 1e-7);
      CHECK(std::abs(norm(p) - R) <= 1e-7);
    }
  }
}
