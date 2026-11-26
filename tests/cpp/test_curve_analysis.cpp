#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <numbers>

#include "nurbforge/analysis/curve_analysis.hpp"
#include "nurbforge/core/error.hpp"
#include "nurbforge/geometry/primitives.hpp"
#include "nurbforge/geometry/transform.hpp"
#include "nurbforge/modeling/blend.hpp"
#include "nurbforge/modeling/surface_tools.hpp"
#include "nurbforge/nurbs/join.hpp"

using namespace nf;
constexpr double pi = std::numbers::pi;

TEST_CASE("curvature and osculating circle", "[curve_analysis]") {
  const Vec3 centre{1, 2, 3};
  const NurbsCurve circle = make_circle(centre, {0, 0, 1}, 2.5);
  for (const double u : {0.0, 0.13, 0.5, 0.77}) {
    const CurveCurvature k = curve_curvature(circle, u);
    CHECK(std::abs(k.curvature - 0.4) < 1e-12);
    CHECK(std::abs(k.radius - 2.5) < 1e-11);
    CHECK(distance(k.center, centre) < 1e-11);
    CHECK(distance(k.point + k.radius * k.normal, centre) < 1e-11);
    CHECK(std::abs(dot(k.tangent, k.normal)) < 1e-12);
    CHECK(distance(k.vector, 0.4 * k.normal) < 1e-12);
  }
  const CurveCurvature flat = curve_curvature(make_line({0, 0, 0}, {3, 1, 0}), 0.4);
  CHECK(flat.curvature == 0.0);
  CHECK(std::isinf(flat.radius));
  CHECK(norm(flat.normal) == 0.0);
  CHECK(distance(flat.center, flat.point) == 0.0);
}

TEST_CASE("curvature graph", "[curve_analysis]") {
  const NurbsCurve circle = make_circle({0, 0, 0}, {0, 0, 1}, 2.0);
  const CurvatureGraph g = curvature_graph(circle, 3.0, 8);
  REQUIRE(g.points.size() == g.tips.size());
  REQUIRE(g.points.size() == g.parameters.size());
  CHECK(std::abs(g.max_curvature - 0.5) < 1e-12);
  for (std::size_t i = 0; i < g.points.size(); ++i) {
    CHECK(std::abs(norm(g.points[i]) - 2.0) < 1e-12);
    CHECK(std::abs(norm(g.tips[i]) - 3.5) < 1e-11);  // outwards by scale * curvature
  }
  // a line running tangentially into an arc: the graph steps at the joint
  const NurbsCurve track = join(make_line({-2, 1, 0}, {0, 1, 0}),
                                make_arc({0, 0, 0}, {0, 1, 0}, {1, 0, 0}, 1.0, 0.0, pi / 2));
  const CurvatureGraph h = curvature_graph(track, 1.0, 4);
  bool step = false;
  for (std::size_t i = 0; i + 1 < h.parameters.size(); ++i)
    if (h.parameters[i] == h.parameters[i + 1] && h.curvatures[i] == 0.0 &&
        std::abs(h.curvatures[i + 1] - 1.0) < 1e-12)
      step = true;
  CHECK(step);
  CHECK_THROWS_AS(curvature_graph(circle, 1.0, 0), GeometryError);
}

TEST_CASE("direction and sense of rotation", "[curve_analysis]") {
  const NurbsCurve ccw = make_circle({0, 0, 0}, {0, 0, 1}, 2.0);
  const CurveDirection d = curve_direction(ccw);
  CHECK(d.closed);
  CHECK(d.planar);
  CHECK(distance(d.plane_normal, {0, 0, 1}) < 1e-12);
  CHECK(d.orientation == 1);
  CHECK(std::abs(d.signed_area - 4 * pi) < 1e-9);
  const CurveDirection r = curve_direction(ccw.reversed());
  CHECK(r.orientation == -1);
  CHECK(distance(r.start_tangent, -1.0 * d.end_tangent) < 1e-12);
  CHECK(curve_direction(ccw, {0, 0, -1}).orientation == -1);  // seen from below
  // a tilted ellipse relative to a given normal
  const NurbsCurve tilted = ccw.transformed(Transform::rotation({0, 0, 0}, {1, 0, 0}, 0.5));
  CHECK(std::abs(curve_direction(tilted, {0, 0, 1}).signed_area - 4 * pi * std::cos(0.5)) < 1e-9);

  const NurbsCurve open = make_arc({0, 0, 0}, {1, 0, 0}, {0, 1, 0}, 1.0, 0.0, pi / 2);
  const CurveDirection o = curve_direction(open);
  CHECK(!o.closed);
  CHECK(o.orientation == 0);
  CHECK(distance(o.start, {1, 0, 0}) < 1e-12);
  CHECK(distance(o.start_tangent, {0, 1, 0}) < 1e-12);
  CHECK(distance(o.end_tangent, {-1, 0, 0}) < 1e-12);
  CHECK(is_closed(ccw));
  CHECK(!is_closed(open));
  CHECK(is_planar(open));
  CHECK(!is_planar(make_helix({0, 0, 0}, {0, 0, 1}, {1, 0, 0}, 1.0, 2.0)));
}

TEST_CASE("continuity between curves and at knots", "[curve_analysis]") {
  const NurbsCurve line = make_line({-2, 1, 0}, {0, 1, 0});
  const NurbsCurve arc = make_arc({0, 0, 0}, {0, 1, 0}, {1, 0, 0}, 1.0, 0.0, pi / 2);
  // tangent, but the curvature jumps from 0 to 1
  const CurveContinuity g1 = continuity(line, CurveEnd::End, arc, CurveEnd::Start);
  CHECK(g1.geometric == 1);
  CHECK(g1.gap < 1e-12);
  CHECK(g1.tangent_angle < 1e-12);
  CHECK(g1.curvature_a == 0.0);
  CHECK(std::abs(g1.curvature_b - 1.0) < 1e-12);
  CHECK(std::abs(g1.curvature_difference - 1.0) < 1e-12);
  // the same joint found from the nearest ends, also with a reversed curve
  CHECK(continuity(line, arc).geometric == 1);
  CHECK(continuity(arc.reversed(), line.reversed()).geometric == 1);
  CHECK(continuity(line.reversed(), CurveEnd::Start, arc, CurveEnd::Start).geometric == 1);
  // a corner, and no contact
  const CurveContinuity g0 = continuity(line, make_line({0, 1, 0}, {0, 3, 0}));
  CHECK(g0.geometric == 0);
  CHECK(std::abs(g0.tangent_angle - pi / 2) < 1e-12);
  CHECK(continuity(line, make_line({0, 1.5, 0}, {0, 3, 0})).geometric == -1);
  // two halves of one curve are C2, a G2 blend joins with curvature continuity
  const NurbsCurve spline =
      NurbsCurve::from_points(3, {{0, 0, 0}, {1, 2, 0}, {3, 2, 1}, {4, 0, 1}, {6, 1, 0}});
  const auto [left, right] = spline.split(0.4);
  const CurveContinuity c2 = continuity(left, right);
  CHECK(c2.geometric == 2);
  CHECK(c2.parametric == 2);

  // knots of one curve
  CHECK(continuity_breaks(spline).empty());
  CHECK(continuity_breaks(make_circle({0, 0, 0}, {0, 0, 1}, 1.0)).empty());
  const NurbsCurve track = join(line, arc);
  const auto breaks = continuity_breaks(track);
  REQUIRE(breaks.size() == 1);
  CHECK(breaks[0].continuity.geometric == 1);
  CHECK(distance(track.evaluate(breaks[0].parameter), {0, 1, 0}) < 1e-12);
  // a rectangle outline: four corners, one of them at the seam
  const auto corners = continuity_breaks(make_rectangle_curve({0, 0, 0}, {2, 0, 0}, {0, 1, 0}));
  REQUIRE(corners.size() == 4);
  for (const auto& b : corners) {
    CHECK(b.continuity.geometric == 0);
    CHECK(std::abs(b.continuity.tangent_angle - pi / 2) < 1e-12);
  }
}

TEST_CASE("duplicate curves", "[curve_analysis]") {
  const NurbsCurve circle = make_circle({0, 0, 0}, {0, 0, 1}, 1.0);
  const NurbsCurve line = make_line({0, 0, 0}, {2, 1, 0});
  const NurbsCurve spline =
      NurbsCurve::from_points(3, {{0, 0, 0}, {1, 2, 0}, {3, 2, 1}, {4, 0, 1}, {6, 1, 0}});
  const RebuiltCurve refit = rebuild_curve(circle, 60);
  REQUIRE(refit.max_error > 1e-9);
  REQUIRE(refit.max_error < 3e-4);
  const std::vector<NurbsCurve> curves{
      circle,                                    // 0
      line,                                      // 1
      spline,                                    // 2
      line.reversed(),                           // 3: the same line the other way
      spline.elevate_degree().insert_knot(0.3),  // 4: the same shape, other representation
      make_circle({0, 0, 0}, {0, 0, 1}, 1.001),  // 5: slightly larger
      refit.curve,                               // 6: a refit of the circle
      make_line({0, 0, 0}, {1, 0.5, 0}),         // 7: half of the line
  };
  const auto tight = find_duplicate_curves(curves, 1e-9);
  REQUIRE(tight.size() == 2);
  CHECK(tight[0] == std::vector<std::size_t>{1, 3});
  CHECK(tight[1] == std::vector<std::size_t>{2, 4});
  const auto loose = find_duplicate_curves(curves, 3 * refit.max_error);
  REQUIRE(loose.size() == 3);
  CHECK(loose[0] == std::vector<std::size_t>{0, 6});
  CHECK(find_duplicate_curves({circle, line}).empty());
  CHECK_THROWS_AS(find_duplicate_curves(curves, 0.0), GeometryError);
}
