#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <numbers>

#include "nurbforge/analysis/properties.hpp"
#include "nurbforge/analysis/silhouette.hpp"
#include "nurbforge/core/error.hpp"
#include "nurbforge/geometry/primitives.hpp"
#include "nurbforge/modeling/construct.hpp"
#include "nurbforge/topology/solids.hpp"

using namespace nf;
constexpr double pi = std::numbers::pi;

namespace {

template <class F>
void for_points(const SilhouetteCurve& c, F f) {
  const double a = c.curve.domain_start(), b = c.curve.domain_end();
  for (int i = 0; i <= 50; ++i) f(a + (b - a) * i / 50.0);
}

double total_length(const std::vector<SilhouetteCurve>& cs) {
  double l = 0.0;
  for (const auto& c : cs) l += arc_length(c.curve);
  return l;
}

// the 3D curve is the image of the p-curve
void check_on_surface(const std::vector<SilhouetteCurve>& cs, const NurbsSurface& s, double tol) {
  for (const auto& c : cs)
    for_points(c, [&](double t) {
      const Vec3 q = c.pcurve.evaluate(t);
      CHECK(distance(s.evaluate(q.x, q.y), c.curve.evaluate(t)) < tol);
    });
}

}  // namespace

TEST_CASE("silhouette of a sphere is a great circle (parallel view)", "[silhouette]") {
  const double R = 2.0;
  const NurbsSurface s = make_sphere({0, 0, 0}, R).faces()[0].surface();
  for (const Vec3& dir : {Vec3{0, 0, -1}, Vec3{1, 0, 0}, Vec3{1, 2, 3}}) {
    const Vec3 d = normalized(dir);
    const auto cs = silhouette(s, View::parallel(dir));
    REQUIRE(!cs.empty());
    for (const auto& c : cs) {
      CHECK(!c.crease);
      CHECK(c.cusps.empty());
      for_points(c, [&](double t) {
        const Vec3 p = c.curve.evaluate(t);
        CHECK(std::abs(norm(p) - R) < 1e-6);
        CHECK(std::abs(dot(p, d)) < 1e-6);
      });
    }
    CHECK(std::abs(total_length(cs) - 2 * pi * R) < 1e-5);
    check_on_surface(cs, s, 2e-6);
  }
}

TEST_CASE("silhouette of a sphere from an eye point is a small circle", "[silhouette]") {
  const double R = 2.0, D = 5.0;
  const NurbsSurface s = make_sphere({0, 0, 0}, R).faces()[0].surface();
  const Vec3 eye = D * normalized(Vec3{1, -2, 0.5});
  const auto cs = silhouette(s, View::perspective(eye));
  REQUIRE(!cs.empty());
  for (const auto& c : cs)
    for_points(c, [&](double t) {
      const Vec3 p = c.curve.evaluate(t);
      CHECK(std::abs(norm(p) - R) < 1e-6);
      CHECK(std::abs(dot(p, eye) - R * R) < 1e-5);  // tangent cone from the eye
    });
  CHECK(std::abs(total_length(cs) - 2 * pi * R * std::sqrt(1 - R * R / (D * D))) < 1e-5);
}

TEST_CASE("silhouettes of a cylinder: straight lines from the side, none along the axis",
          "[silhouette]") {
  const double R = 1.5, H = 4.0;
  const NurbsSurface s = extrude(make_circle({0, 0, 0}, {0, 0, 1}, R), {0, 0, H});
  const auto cs = silhouette(s, View::parallel({1, 0, 0}));
  REQUIRE(cs.size() == 2);
  for (const auto& c : cs) {
    for_points(c, [&](double t) {
      const Vec3 p = c.curve.evaluate(t);
      CHECK(std::abs(p.x) < 1e-9);
      CHECK(std::abs(std::abs(p.y) - R) < 1e-9);
    });
    CHECK(std::abs(arc_length(c.curve) - H) < 1e-9);
  }
  CHECK(cs[0].curve.evaluate(cs[0].curve.domain_start()).y *
            cs[1].curve.evaluate(cs[1].curve.domain_start()).y <
        0.0);
  CHECK(silhouette(s, View::parallel({0, 0, 1})).empty());

  // clipped to a trimmed face: the lower half of the cylinder
  const auto du = s.domain_u(), dv = s.domain_v();
  const TrimmedSurface half(s, rectangle_loop(du[0], du[1], dv[0], 0.5 * (dv[0] + dv[1])));
  const auto hs = silhouette(half, View::parallel({1, 1, 0}));  // away from the seam
  REQUIRE(hs.size() == 2);
  CHECK(std::abs(total_length(hs) - H) < 1e-7);

  // a box has no silhouettes, a solid cylinder only those of its side
  CHECK(silhouette(make_box({0, 0, 0}, {1, 2, 3}), View::parallel({1, 2, 3})).empty());
  const Shell solid = make_cylinder({0, 0, 0}, {0, 0, 1}, R, H);
  const auto ss = silhouette(solid, View::parallel({1, 1, 0}));
  CHECK(ss.size() == 2);
  CHECK(std::abs(total_length(ss) - 2 * H) < 1e-7);
  // silhouettes on the seam of the side face and on a knot line
  for (const Vec3& dir : {Vec3{0, 1, 0}, Vec3{1, 0, 0}, Vec3{0, -1, 0}, Vec3{-1, 0, 0}}) {
    const auto axis =
        silhouette(make_cylinder({3, 1, 1}, {0, 0, 1}, 0.4, 1.5), View::parallel(dir));
    REQUIRE(axis.size() == 2);
    CHECK(std::abs(total_length(axis) - 3.0) < 1e-9);
  }
}

TEST_CASE("silhouettes of a torus: circles along the axis, cusps in an oblique view",
          "[silhouette]") {
  const double R = 3.0, r = 1.0;
  const NurbsSurface s = make_torus({0, 0, 0}, {0, 0, 1}, R, r).faces()[0].surface();
  const auto top = silhouette(s, View::parallel({0, 0, -1}));
  REQUIRE(!top.empty());
  for (const auto& c : top)
    for_points(c, [&](double t) {
      const Vec3 p = c.curve.evaluate(t);
      CHECK(std::abs(p.z) < 1e-6);
      CHECK(std::abs(std::abs(std::hypot(p.x, p.y) - R) - r) < 1e-6);
    });
  CHECK(std::abs(total_length(top) - 4 * pi * R) < 1e-4);

  // oblique: every point has a normal perpendicular to the view; the inner fold has four cusps
  const Vec3 d = normalized(Vec3{std::sin(1.0), 0, -std::cos(1.0)});
  const auto cs = silhouette(s, View::parallel(d));
  std::size_t cusps = 0;
  for (const auto& c : cs) {
    cusps += c.cusps.size();
    for_points(c, [&](double t) {
      const Vec3 q = c.pcurve.evaluate(t);
      CHECK(std::abs(dot(s.normal(q.x, q.y), d)) < 1e-5);
    });
  }
  CHECK(cusps == 4);
  check_on_surface(cs, s, 2e-6);
}

TEST_CASE("a crease between a facing and an averted part is a silhouette", "[silhouette]") {
  // roof: two planar parts meeting in the ridge x = 0, z = 1
  const NurbsSurface roof = extrude(make_polyline({{-1, 0, 0}, {0, 0, 1}, {1, 0, 0}}), {0, 3, 0});
  const auto side = silhouette(roof, View::parallel({1, 0, 0}));
  REQUIRE(side.size() == 1);
  CHECK(side[0].crease);
  for_points(side[0], [&](double t) {
    const Vec3 p = side[0].curve.evaluate(t);
    CHECK(std::abs(p.x) < 1e-12);
    CHECK(std::abs(p.z - 1.0) < 1e-12);
  });
  CHECK(std::abs(arc_length(side[0].curve) - 3.0) < 1e-9);
  CHECK(silhouette(roof, View::parallel({0, 0, -1})).empty());  // both parts face the viewer
  CHECK_THROWS_AS(silhouette(roof, View::parallel({0, 0, 0})), GeometryError);
}
