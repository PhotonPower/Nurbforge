#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <numbers>
#include <vector>

#include "nurbforge/core/error.hpp"
#include "nurbforge/geometry/primitives.hpp"
#include "nurbforge/modeling/construct.hpp"
#include "nurbforge/modeling/gordon.hpp"
#include "nurbforge/modeling/surface_tools.hpp"
#include "nurbforge/modeling/sweep.hpp"
#include "nurbforge/nurbs/knot_vector.hpp"
#include "nurbforge/topology/shell.hpp"
#include "nurbforge/topology/solid_ops.hpp"
#include "nurbforge/topology/trimmed_surface.hpp"

using namespace nf;
constexpr double pi = std::numbers::pi;

namespace {
double radius_xy(const Vec3& p) { return std::hypot(p.x, p.y); }
}  // namespace

TEST_CASE("rebuild curves and surfaces", "[surface_tools]") {
  SECTION("a spline from the target space is reproduced") {
    const NurbsCurve c(
        3, clamped_uniform_knots(3, 8),
        {{0, 0, 0}, {1, 2, 0}, {2, -1, 1}, {3, 0, 2}, {4, 3, 0}, {5, 1, 1}, {6, 0, 0}, {7, 2, 2}});
    const RebuiltCurve r = rebuild_curve(c, 8, 3);
    CHECK(r.max_error <= 1e-12);
    CHECK(r.curve.num_control_points() == 8);
  }
  SECTION("a circle: error reported and shrinking with more control points") {
    const NurbsCurve circle = make_circle({0, 0, 0}, {0, 0, 1}, 1.0);
    const RebuiltCurve a = rebuild_curve(circle, 12), b = rebuild_curve(circle, 40);
    CHECK(!a.curve.is_rational());
    CHECK(a.max_error < 1e-2);
    CHECK(b.max_error < a.max_error / 20);
    // the reported error is the real one
    double e = 0;
    for (int i = 0; i <= 1000; ++i) {
      const double t = i / 1000.0;
      e = std::max(e, distance(circle.evaluate(t), b.curve.evaluate(t)));
    }
    CHECK(e <= b.max_error * (1 + 1e-6));
    CHECK_THROWS_AS(rebuild_curve(circle, 3, 3), GeometryError);
  }
  SECTION("surfaces") {
    const NurbsSurface patch = make_bilinear({0, 0, 0}, {2, 0, 1}, {0, 1, 0}, {2, 1, 3});
    const RebuiltSurface r = rebuild_surface(patch, 5, 6, 3, 2);
    CHECK(r.max_error <= 1e-12);
    CHECK(r.surface.num_u() == 5);
    CHECK(r.surface.num_v() == 6);
    const NurbsSurface sph = revolve(make_arc({0, 0, 0}, {0, 0, -1}, {1, 0, 0}, 1.0, 0.3, pi - 0.3),
                                     {0, 0, 0}, {0, 0, 1}, pi);
    const RebuiltSurface s = rebuild_surface(sph, 12, 12);
    CHECK(s.max_error < 1e-3);
  }
}

TEST_CASE("tapered extrusion, ribbon, extrusion to a point", "[surface_tools]") {
  SECTION("square profile: four planar faces, frustum volume") {
    const NurbsCurve sq = make_rectangle_curve({0, 0, 0}, {2, 0, 0}, {0, 2, 0});
    const double h = 1.0, ang = 10 * pi / 180, t = std::tan(ang);
    const auto faces = extrude_tapered(sq, {0, 0, h}, ang);
    REQUIRE(faces.size() == 4);
    std::vector<TrimmedSurface> tf;
    for (const auto& f : faces) tf.emplace_back(f);
    const Shell s = cap(Shell::sew(tf));
    REQUIRE(s.is_closed());
    const double a1 = 4, a2 = (2 + 2 * t) * (2 + 2 * t);
    CHECK(std::abs(std::abs(s.volume()) - h / 3 * (a1 + a2 + std::sqrt(a1 * a2))) <= 1e-9);
    // negative angle: inwards
    const auto in = extrude_tapered(sq, {0, 0, h}, -ang);
    std::vector<TrimmedSurface> ti;
    for (const auto& f : in) ti.emplace_back(f);
    const double b2 = (2 - 2 * t) * (2 - 2 * t);
    CHECK(std::abs(std::abs(cap(Shell::sew(ti)).volume()) -
                   h / 3 * (a1 + b2 + std::sqrt(a1 * b2))) <= 1e-9);
  }
  SECTION("circle: one conical surface") {
    const NurbsCurve c = make_circle({0, 0, 0}, {0, 0, 1}, 1.0);
    const auto f = extrude_tapered(c, {0, 0, 2}, pi / 12);
    REQUIRE(f.size() == 1);
    for (int i = 0; i <= 20; ++i) {
      const double u = i / 20.0;
      CHECK(std::abs(radius_xy(f[0].evaluate(u, 1.0)) - (1 + 2 * std::tan(pi / 12))) <= 1e-6);
      CHECK(std::abs(f[0].evaluate(u, 1.0).z - 2) <= 1e-9);
    }
    CHECK_THROWS_AS(extrude_tapered(c, {0, 0, 1}, pi / 2), GeometryError);
  }
  SECTION("ribbon along an L: mitred trapezoids") {
    const NurbsCurve l = make_polyline({{0, 0, 0}, {2, 0, 0}, {2, 1, 0}});
    const auto r = ribbon(l, 0.5, {0, 0, 1});
    REQUIRE(r.size() == 2);
    double area = 0;
    for (const auto& s : r) area += TrimmedSurface(s).area();
    CHECK(std::abs(area - 1.25) <= 1e-9);
  }
  SECTION("extrusion to a point: exact cone") {
    const NurbsCurve c = make_circle({0, 0, 0}, {0, 0, 1}, 1.0);
    const NurbsSurface cone = extrude_to_point(c, {0, 0, 2});
    for (int i = 0; i <= 16; ++i)
      CHECK(std::abs(radius_xy(cone.evaluate(i / 16.0, 0.25)) - 0.75) <= 1e-12);
    CHECK(std::abs(TrimmedSurface(cone).area() - pi * std::sqrt(5.0)) <= 1e-8);
  }
}

TEST_CASE("Coons patches with rational boundaries and triangles", "[surface_tools]") {
  SECTION("exact cylinder patch from two arcs and two lines") {
    const NurbsCurve a0 = make_arc({0, 0, 0}, {1, 0, 0}, {0, 1, 0}, 1.0, 0, pi / 2);
    const NurbsCurve a1 = make_arc({0, 0, 1}, {1, 0, 0}, {0, 1, 0}, 1.0, 0, pi / 2);
    const NurbsSurface s =
        coons(a0, a1, make_line({1, 0, 0}, {1, 0, 1}), make_line({0, 1, 0}, {0, 1, 1}));
    CHECK(!s.weights().empty());
    for (int i = 0; i <= 10; ++i)
      for (int j = 0; j <= 10; ++j)
        CHECK(std::abs(radius_xy(s.evaluate(i / 10.0, j / 10.0)) - 1) <= 1e-12);
  }
  SECTION("triangle patch from three arcs of the unit sphere") {
    const NurbsCurve c0 = make_arc({0, 0, 0}, {1, 0, 0}, {0, 1, 0}, 1.0, 0, pi / 2);
    const NurbsCurve c1 = make_arc({0, 0, 0}, {0, 1, 0}, {0, 0, 1}, 1.0, 0, pi / 2);
    const NurbsCurve c2 = make_arc({0, 0, 0}, {0, 0, 1}, {1, 0, 0}, 1.0, 0, pi / 2);
    const NurbsSurface t = triangle_patch(c0, c2, c1);  // any order of the other two
    // the three sides are the exact arcs, the opposite corner of c0 is the pole (0, 0, 1)
    for (int i = 0; i <= 10; ++i) {
      const double u = i / 10.0;
      CHECK(std::abs(norm(t.evaluate(u, 0.0)) - 1) <= 1e-12);
      CHECK(std::abs(norm(t.evaluate(0.0, u)) - 1) <= 1e-12);
      CHECK(std::abs(norm(t.evaluate(1.0, u)) - 1) <= 1e-12);
      CHECK(distance(t.evaluate(u, 1.0), {0, 0, 1}) <= 1e-12);
    }
    const Vec3 n = t.normal(0.5, 1.0);  // limit normal at the degenerate side
    CHECK(std::abs(norm(n) - 1) <= 1e-9);
    CHECK_THROWS_AS(triangle_patch(c0, c1, make_line({0, 0, 2}, {1, 0, 0})), GeometryError);
  }
}

TEST_CASE("sweep with several profiles and rail revolve", "[surface_tools]") {
  SECTION("two circles along a line: frustum") {
    const NurbsCurve path = make_line({0, 0, 0}, {0, 0, 4});
    const std::vector<NurbsCurve> prof{make_circle({0, 0, 0}, {0, 0, 1}, 1.0),
                                       make_circle({0, 0, 4}, {0, 0, 1}, 2.0)};
    SweepOptions o;
    o.tolerance = 1e-8;
    const NurbsSurface s = sweep(prof, path, {}, o);
    for (int i = 0; i <= 8; ++i)
      for (const double v : {0.0, 1.0, 2.0, 3.0, 4.0}) {
        const Vec3 p = s.evaluate(i / 8.0, v / 4);  // the line's parameter is z / 4
        CHECK(std::abs(radius_xy(p) - (1 + v / 4)) <= 1e-7);
        CHECK(std::abs(p.z - v) <= 1e-7);
      }
    // three profiles: the middle one is passed through
    const std::vector<NurbsCurve> three{prof[0], make_circle({0, 0, 2}, {0, 0, 1}, 0.5), prof[1]};
    const NurbsSurface w = sweep(three, path, {0.0, 0.5, 1.0}, o);
    CHECK(std::abs(radius_xy(w.evaluate(0.3, 0.5)) - 0.5) <= 1e-7);
    CHECK_THROWS_AS(sweep(three, path, {0.0, 0.0, 1.0}, o), GeometryError);
  }
  SECTION("rail revolve: a disc with an elliptic rim") {
    const NurbsCurve prof = make_line({0, 0, 0}, {2, 0, 0});
    const NurbsCurve rail = make_circle({0, 0, 0}, {0, 0, 1}, 1.0)
                                .transformed(Transform::scaling({0, 0, 0}, Vec3{2, 1, 1}));
    SweepOptions o;
    o.tolerance = 1e-8;
    const NurbsSurface s = rail_revolve(prof, rail, {0, 0, 0}, {0, 0, 1}, false, o);
    for (int j = 0; j <= 16; ++j) {
      const double v = j / 16.0;
      CHECK(distance(s.evaluate(1.0, v), rail.evaluate(v)) <= 1e-8);
      CHECK(norm(s.evaluate(0.0, v)) <= 1e-12);
    }
    CHECK(std::abs(TrimmedSurface(s).area() - 2 * pi) <= 1e-6);
    // scale_height: a cone whose rim height follows the rail
    const NurbsCurve tilted =
        rail.transformed(Transform({1, 0, 0, 0, 0, 1, 0, 0, 0.25, 0, 0, 1}));  // z = 1 + x/4
    const NurbsSurface c =
        rail_revolve(make_line({0, 0, 0}, {2, 0, 1.5}), tilted, {0, 0, 0}, {0, 0, 1}, true, o);
    for (int j = 0; j <= 16; ++j)
      CHECK(distance(c.evaluate(1.0, j / 16.0), tilted.evaluate(j / 16.0)) <= 1e-8);
  }
}
