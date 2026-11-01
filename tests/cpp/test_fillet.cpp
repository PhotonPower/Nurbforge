#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <algorithm>
#include <cmath>
#include <numbers>

#include "nurbforge/core/error.hpp"
#include "nurbforge/geometry/primitives.hpp"
#include "nurbforge/modeling/construct.hpp"
#include "nurbforge/modeling/fillet.hpp"
#include "nurbforge/nurbs/join.hpp"

using namespace nf;
constexpr double pi = std::numbers::pi;

namespace {

constexpr Vec3 ez{0, 0, 1};

double max_radius_error(const NurbsCurve& arc, const Vec3& centre, double r) {
  double e = 0.0;
  for (int i = 0; i <= 50; ++i) {
    const double t = arc.domain_start() + (arc.domain_end() - arc.domain_start()) * i / 50.0;
    e = std::max(e, std::abs(distance(arc.evaluate(t), centre) - r));
  }
  return e;
}

double angle(const Vec3& a, const Vec3& b) { return std::atan2(norm(cross(a, b)), dot(a, b)); }

// sin of the angle between two (unoriented) directions
double sin_angle(const Vec3& a, const Vec3& b) { return norm(cross(normalized(a), normalized(b))); }

template <class F>
void for_grid(const NurbsSurface& s, int n, F f) {
  const auto du = s.domain_u(), dv = s.domain_v();
  for (int i = 0; i <= n; ++i)
    for (int j = 0; j <= n; ++j)
      f(du[0] + (du[1] - du[0]) * i / n, dv[0] + (dv[1] - dv[0]) * j / n);
}

// Common checks of a rolling-ball fillet: every surface point at distance r from the spine
// point of its section, the surface tangent to s1 / s2 along the contact curves, the contact
// curves equal to the mapped p-curves and the section planes normal to the spine.
void check_rolling_ball(const SurfaceFillet& f, const NurbsSurface& s1, const NurbsSurface& s2,
                        double r, double tol) {
  for_grid(f.surface, 24, [&](double u, double v) {
    REQUIRE(std::abs(distance(f.surface.evaluate(u, v), f.spine.evaluate(v)) - r) <= tol);
  });
  const auto dv = f.surface.domain_v();
  for (int j = 0; j <= 24; ++j) {
    const double t = dv[0] + (dv[1] - dv[0]) * j / 24.0;
    const Vec3 q1 = f.pcurve1.evaluate(t), q2 = f.pcurve2.evaluate(t);
    REQUIRE(distance(s1.evaluate(q1.x, q1.y), f.contact1.evaluate(t)) <= 2 * tol);
    REQUIRE(distance(s2.evaluate(q2.x, q2.y), f.contact2.evaluate(t)) <= 2 * tol);
    REQUIRE(distance(f.contact1.evaluate(t), f.surface.evaluate(0.0, t)) <= 1e-12);
    REQUIRE(distance(f.contact2.evaluate(t), f.surface.evaluate(1.0, t)) <= 1e-12);
    // tangency: fillet normal parallel to the surface normals along the contact curves
    REQUIRE(sin_angle(f.surface.normal(0.0, t), s1.normal(q1.x, q1.y)) <= 1e-4);
    REQUIRE(sin_angle(f.surface.normal(1.0, t), s2.normal(q2.x, q2.y)) <= 1e-4);
    // the section arc lies in the plane normal to the spine tangent
    const Vec3 T = normalized(f.spine.derivatives(t, 1)[1]);
    REQUIRE(std::abs(dot(f.surface.evaluate(0.5, t) - f.spine.evaluate(t), T)) <= 10 * tol);
  }
}

NurbsSurface sphere(const Vec3& c, double R) {
  return revolve(make_arc(c, {1, 0, 0}, {0, 0, 1}, R, -pi / 2, pi / 2), c, ez, 2 * pi);
}

// +1 if the surface normal at (u, v) points away from `axis_point` (orthogonally to `axis`)
int outward(const NurbsSurface& s, double u, double v, const Vec3& axis_point, const Vec3& axis) {
  Vec3 d = s.evaluate(u, v) - axis_point;
  d = d - dot(d, axis) * axis;
  return dot(s.normal(u, v), d) > 0.0 ? 1 : -1;
}

}  // namespace

TEST_CASE("fillet between two lines", "[fillet]") {
  const double r = 0.4;
  const double turn = GENERATE(30.0, 90.0, 135.0, -60.0, -120.0);  // degrees, left > 0
  const double phi = turn * pi / 180.0;
  const Vec3 X{1, 2, 0};
  const auto a = make_line({-2, 2, 0}, X);
  const auto b = make_line(X, X + 3.0 * Vec3{std::cos(phi), std::sin(phi), 0});
  const auto f = fillet_curves(a, b, r, ez);
  const double theta = pi - std::abs(phi);  // angle between the kept parts
  const double setback = r / std::tan(theta / 2.0);
  CHECK(std::abs(distance(f.point_a, X) - setback) <= 1e-12);
  CHECK(std::abs(distance(f.point_b, X) - setback) <= 1e-12);
  CHECK(f.point_a.x < X.x);
  CHECK(std::abs(distance(f.center, X) - r / std::sin(theta / 2.0)) <= 1e-12);
  // centre on the inner side of the corner, at distance r from both lines
  CHECK(std::abs(f.center.y - X.y - (phi > 0 ? r : -r)) <= 1e-12);
  // exact rational arc
  CHECK(f.arc.degree() == 2);
  CHECK(f.arc.is_rational());
  CHECK(max_radius_error(f.arc, f.center, r) <= 1e-12);
  CHECK(distance(f.arc.evaluate(0.0), f.point_a) <= 1e-12);
  CHECK(distance(f.arc.evaluate(1.0), f.point_b) <= 1e-12);
  CHECK(std::abs(angle(f.point_a - f.center, f.point_b - f.center) - std::abs(phi)) <= 1e-12);
  // a, arc, b form a G1 chain
  REQUIRE(f.trimmed_a);
  REQUIRE(f.trimmed_b);
  CHECK(distance(f.trimmed_a->evaluate(f.trimmed_a->domain_start()), a.evaluate(0.0)) <= 1e-12);
  CHECK(distance(f.trimmed_b->evaluate(f.trimmed_b->domain_end()), b.evaluate(1.0)) <= 1e-12);
  CHECK(angle(f.trimmed_a->tangent(f.trimmed_a->domain_end()), f.arc.tangent(0.0)) <= 1e-9);
  CHECK(angle(f.arc.tangent(1.0), f.trimmed_b->tangent(f.trimmed_b->domain_start())) <= 1e-9);
  const auto chain = join(join(*f.trimmed_a, f.arc, 1e-9), *f.trimmed_b, 1e-9);
  CHECK(distance(chain.evaluate(chain.domain_start()), a.evaluate(0.0)) <= 1e-9);
  CHECK(distance(chain.evaluate(chain.domain_end()), b.evaluate(1.0)) <= 1e-9);
}

TEST_CASE("fillet of reversed lines keeps the corner", "[fillet]") {
  // a ends at the corner, b also ends there: kept parts are both before the corner
  const auto a = make_line({0, 0, 0}, {2, 0, 0});
  const auto b = make_line({2, 2, 0}, {2, 0, 0});
  const auto f = fillet_curves(a, b, 0.5, ez);
  CHECK(distance(f.center, {1.5, 0.5, 0}) <= 1e-12);
  REQUIRE(f.trimmed_b);
  CHECK(distance(f.trimmed_b->evaluate(f.trimmed_b->domain_start()), {2, 0.5, 0}) <= 1e-12);
  CHECK(distance(f.trimmed_b->evaluate(f.trimmed_b->domain_end()), {2, 2, 0}) <= 1e-12);
}

TEST_CASE("fillet between a line and a circle is analytic", "[fillet]") {
  const double R = 2.0, h = 0.5, r = 0.3;
  const auto a = make_line({-4, h, 0}, {4, h, 0});
  const auto b = make_arc({0, 0, 0}, {1, 0, 0}, {0, 1, 0}, R, -pi / 2, pi / 2);

  SECTION("automatic side: corner between the kept parts (inside the circle, below the line)") {
    const auto f = fillet_curves(a, b, r, ez);
    const double cx = std::sqrt((R - r) * (R - r) - (h - r) * (h - r));
    CHECK(distance(f.center, {cx, h - r, 0}) <= 1e-10);
    CHECK(distance(f.point_a, {cx, h, 0}) <= 1e-10);
    CHECK(distance(f.point_b, R / (R - r) * f.center) <= 1e-10);
    CHECK(max_radius_error(f.arc, f.center, r) <= 1e-12);
    REQUIRE(f.trimmed_a);
    REQUIRE(f.trimmed_b);
    CHECK(angle(f.trimmed_a->tangent(f.trimmed_a->domain_end()), f.arc.tangent(0.0)) <= 1e-8);
    CHECK(angle(f.arc.tangent(1.0), f.trimmed_b->tangent(f.trimmed_b->domain_start())) <= 1e-8);
  }
  SECTION("explicit sides: above the line, outside the circle") {
    CurveFilletOptions opt;
    opt.side_a = 1;   // left of +x: above
    opt.side_b = -1;  // right of the counter-clockwise arc: outside
    const auto f = fillet_curves(a, b, r, ez, opt);
    const double cx = std::sqrt((R + r) * (R + r) - (h + r) * (h + r));
    CHECK(distance(f.center, {cx, h + r, 0}) <= 1e-10);
    CHECK(distance(f.point_a, {cx, h, 0}) <= 1e-10);
    CHECK(distance(f.point_b, R / (R + r) * f.center) <= 1e-10);
    CHECK(max_radius_error(f.arc, f.center, r) <= 1e-12);
    REQUIRE(f.trimmed_a);
    CHECK(angle(f.trimmed_a->tangent(f.trimmed_a->domain_end()), f.arc.tangent(0.0)) <= 1e-8);
    CHECK(angle(f.arc.tangent(1.0), f.trimmed_b->tangent(f.trimmed_b->domain_start())) <= 1e-8);
  }
}

TEST_CASE("fillet_curves rejects invalid input", "[fillet]") {
  const auto a = make_line({0, 0, 0}, {1, 0, 0});
  const auto b = make_line({1, 0, 0}, {1, 1, 0});
  CHECK_THROWS_AS(fillet_curves(a, b, 0.0, ez), GeometryError);
  CHECK_THROWS_AS(fillet_curves(a, b, 2.0, ez), GeometryError);  // too large for the segments
  CHECK_THROWS_AS(fillet_curves(a, make_line({0, 1, 0}, {1, 1, 0}), 0.2, ez), GeometryError);
  CHECK_THROWS_AS(fillet_curves(a, b, 0.2, {1, 0, 0}), GeometryError);  // not in that plane
  CurveFilletOptions bad;
  bad.side_a = 2;
  CHECK_THROWS_AS(fillet_curves(a, b, 0.2, ez, bad), GeometryError);
}

TEST_CASE("rolling-ball fillet between perpendicular planes is a quarter cylinder", "[fillet]") {
  const double r = 0.5;
  const auto s1 = make_rectangle({-1, 0, 0}, {5, 0, 0}, {0, 3, 0});   // z = 0, normal +z
  const auto s2 = make_rectangle({0, 0, -1}, {0, 3, 0}, {0, 0, 5});   // x = 0, normal +x
  SurfaceFilletOptions opt;
  opt.tolerance = 1e-7;
  const auto f = fillet_surfaces(s1, s2, r, opt);
  check_rolling_ball(f, s1, s2, r, 1e-7);
  // spine: the line x = z = r over the full width
  for (int j = 0; j <= 20; ++j) {
    const auto dv = f.surface.domain_v();
    const double t = dv[0] + (dv[1] - dv[0]) * j / 20.0;
    const Vec3 c = f.spine.evaluate(t), p1 = f.contact1.evaluate(t), p2 = f.contact2.evaluate(t);
    CHECK(std::abs(c.x - r) <= 1e-9);
    CHECK(std::abs(c.z - r) <= 1e-9);
    CHECK(std::abs(p1.z) <= 1e-9);
    CHECK(std::abs(p1.x - r) <= 1e-9);
    CHECK(std::abs(p2.x) <= 1e-9);
    CHECK(std::abs(p2.z - r) <= 1e-9);
  }
  const Vec3 e0 = f.spine.evaluate(f.spine.domain_start()), e1 = f.spine.evaluate(f.spine.domain_end());
  CHECK(std::abs(std::min(e0.y, e1.y)) <= 1e-9);
  CHECK(std::abs(std::max(e0.y, e1.y) - 3.0) <= 1e-9);
  // the quarter facing the edge: every point within the corner x <= r, z <= r, at distance r
  // from the spine line
  for_grid(f.surface, 20, [&](double u, double v) {
    const Vec3 p = f.surface.evaluate(u, v);
    CHECK(std::abs(std::hypot(p.x - r, p.z - r) - r) <= 1e-7);
    CHECK(p.x <= r + 1e-9);
    CHECK(p.z <= r + 1e-9);
  });
}

TEST_CASE("rolling-ball fillet between a plane and a perpendicular cylinder", "[fillet]") {
  const double R = 1.0, r = 0.4, tol = 1e-6;
  const auto plane = make_rectangle({-3, -3, 0}, {6, 0, 0}, {0, 6, 0});  // normal +z
  const auto cyl = extrude(make_circle({0, 0, -1}, ez, R), {0, 0, 3});
  SurfaceFilletOptions opt;
  opt.tolerance = tol;
  opt.side1 = 1;                                          // above the plane
  opt.side2 = outward(cyl, 0.3, 0.5, {0, 0, 0}, ez);     // outside the cylinder
  const auto f = fillet_surfaces(plane, cyl, r, opt);
  check_rolling_ball(f, plane, cyl, r, tol);
  // the spine is the full circle (from seam to seam)
  CHECK(distance(f.spine.evaluate(f.spine.domain_start()), f.spine.evaluate(f.spine.domain_end())) <= tol);
  const auto dv = f.surface.domain_v();
  for (int j = 0; j <= 40; ++j) {
    const double t = dv[0] + (dv[1] - dv[0]) * j / 40.0;
    const Vec3 c = f.spine.evaluate(t), p1 = f.contact1.evaluate(t), p2 = f.contact2.evaluate(t);
    CHECK(std::abs(std::hypot(c.x, c.y) - (R + r)) <= tol);
    CHECK(std::abs(c.z - r) <= tol);
    CHECK(std::abs(std::hypot(p1.x, p1.y) - (R + r)) <= tol);
    CHECK(std::abs(p1.z) <= tol);
    CHECK(std::abs(std::hypot(p2.x, p2.y) - R) <= tol);
    CHECK(std::abs(p2.z - r) <= tol);
  }
}

TEST_CASE("rolling-ball fillet between a plane and a sphere", "[fillet]") {
  const double R = 1.0, r = 0.3, tol = 1e-6;
  const Vec3 O{0, 0, 0.5};
  const auto plane = make_rectangle({-3, -3, 0}, {6, 0, 0}, {0, 6, 0});  // normal +z
  const auto sph = sphere(O, R);
  SurfaceFilletOptions opt;
  opt.tolerance = tol;
  opt.offset_tolerance = 1e-5;  // start values only
  opt.side2 = dot(sph.normal(0.3, 0.4), sph.evaluate(0.3, 0.4) - O) > 0 ? 1 : -1;  // outside
  const auto f = fillet_surfaces(plane, sph, r, opt);
  check_rolling_ball(f, plane, sph, r, tol);
  CHECK(distance(f.spine.evaluate(f.spine.domain_start()), f.spine.evaluate(f.spine.domain_end())) <= tol);
  const auto dv = f.surface.domain_v();
  for (int j = 0; j <= 40; ++j) {
    const double t = dv[0] + (dv[1] - dv[0]) * j / 40.0;
    const Vec3 c = f.spine.evaluate(t), p2 = f.contact2.evaluate(t);
    CHECK(std::abs(distance(c, O) - (R + r)) <= tol);
    CHECK(std::abs(c.z - r) <= tol);
    CHECK(std::abs(distance(p2, O) - R) <= tol);
    CHECK(std::abs(f.contact1.evaluate(t).z) <= tol);
  }
}

TEST_CASE("fillet_surfaces rejects invalid input", "[fillet]") {
  const auto s1 = make_rectangle({-1, 0, 0}, {5, 0, 0}, {0, 3, 0});
  const auto s2 = make_rectangle({0, 0, -1}, {0, 3, 0}, {0, 0, 5});
  CHECK_THROWS_AS(fillet_surfaces(s1, s2, -1.0), GeometryError);
  SurfaceFilletOptions opt;
  opt.side1 = 0;
  CHECK_THROWS_AS(fillet_surfaces(s1, s2, 0.5, opt), GeometryError);
  // parallel planes: the offsets never meet
  const auto s3 = make_rectangle({-1, 0, 2}, {5, 0, 0}, {0, 3, 0});
  CHECK_THROWS_AS(fillet_surfaces(s1, s3, 0.5), GeometryError);
}
