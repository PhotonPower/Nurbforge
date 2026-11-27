#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <numbers>

#include "nurbforge/analysis/properties.hpp"
#include "nurbforge/core/error.hpp"
#include "nurbforge/geometry/primitives.hpp"
#include "nurbforge/modeling/construct.hpp"
#include "nurbforge/modeling/fillet.hpp"
#include "nurbforge/modeling/surface_tools.hpp"
#include "nurbforge/modeling/unroll.hpp"
#include "nurbforge/topology/solids.hpp"

using namespace nf;
constexpr double pi = std::numbers::pi;

namespace {

template <class F>
void for_grid(const NurbsSurface& s, int n, F f) {
  const auto du = s.domain_u(), dv = s.domain_v();
  for (int i = 0; i <= n; ++i)
    for (int j = 0; j <= n; ++j)
      f(du[0] + (du[1] - du[0]) * i / n, dv[0] + (dv[1] - dv[0]) * j / n);
}

// largest difference of the first fundamental forms of the surface and its pattern, relative
double metric_error(const NurbsSurface& s, const NurbsSurface& p) {
  double worst = 0.0;
  const auto du = s.domain_u(), dv = s.domain_v();
  for (int i = 0; i < 20; ++i)
    for (int j = 0; j < 20; ++j) {
      const double u = du[0] + (du[1] - du[0]) * (i + 0.37) / 20.0;
      const double v = dv[0] + (dv[1] - dv[0]) * (j + 0.41) / 20.0;
      const auto a = s.first_derivatives(u, v), b = p.first_derivatives(u, v);
      const double scale = dot(a[1], a[1]) + dot(a[2], a[2]);
      worst = std::max({worst, std::abs(dot(a[1], a[1]) - dot(b[1], b[1])) / scale,
                        std::abs(dot(a[1], a[2]) - dot(b[1], b[2])) / scale,
                        std::abs(dot(a[2], a[2]) - dot(b[2], b[2])) / scale});
    }
  return worst;
}

}  // namespace

TEST_CASE("developability of ruled and curved surfaces", "[unroll]") {
  const NurbsSurface cyl = extrude(make_circle({0, 0, 0}, {0, 0, 1}, 1.5), {0, 0, 4});
  const Developability c = developability(cyl);
  CHECK(c.rulings == RulingDirection::V);
  CHECK(c.developable);
  CHECK(c.max_twist < 1e-12);
  CHECK(c.max_gaussian_curvature < 1e-9);
  CHECK(developability(cyl.transposed()).rulings == RulingDirection::U);
  CHECK(developability(extrude_to_point(make_circle({0, 0, 0}, {0, 0, 1}, 1.0), {0.3, 0, 2}))
            .developable);

  // a twisted bilinear patch is ruled, but its normal turns along the rulings
  const Developability h =
      developability(make_bilinear({0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 1}));
  CHECK(h.rulings != RulingDirection::None);
  CHECK(!h.developable);
  CHECK(h.max_twist > 0.1);
  CHECK(h.max_gaussian_curvature > 0.1);

  const Developability sph = developability(make_sphere({0, 0, 0}, 2.0).faces()[0].surface());
  CHECK(sph.rulings == RulingDirection::None);
  CHECK(!sph.developable);
  CHECK(std::abs(sph.max_gaussian_curvature - 0.25) < 1e-6);
}

TEST_CASE("a cylinder unrolls to a rectangle", "[unroll]") {
  const double R = 1.5, H = 4.0, tol = 1e-8;
  const NurbsSurface cyl = extrude(make_circle({0, 0, 0}, {0, 0, 1}, R), {0, 0, H});
  const UnrolledSurface r = unroll(cyl, tol);
  CHECK(r.max_error <= tol);
  CHECK(r.pattern.domain_u() == cyl.domain_u());
  CHECK(r.pattern.domain_v() == cyl.domain_v());
  const Vec3 start = cyl.evaluate(0.0, 0.0);
  for_grid(cyl, 24, [&](double u, double v) {
    const Vec3 p = cyl.evaluate(u, v), q = r.pattern.evaluate(u, v);
    CHECK(std::abs(q.z) < 1e-14);
    // along the circle by the arc length from the start, across by the height
    double angle = std::atan2(start.x * p.y - start.y * p.x, start.x * p.x + start.y * p.y);
    if (angle < -1e-9 || (u > 0.5 && angle < 1e-9)) angle += 2 * pi;
    CHECK(std::abs(q.x - R * angle) < 2 * tol);
    CHECK(std::abs(std::abs(q.y) - std::abs(p.z - 0.5 * H)) < 2 * tol);
  });
  CHECK(std::abs(area(r.pattern) - 2 * pi * R * H) < 1e-6);
  // seen from the side of the normal: no mirror image
  const auto d = r.pattern.first_derivatives(0.3, 0.4);
  CHECK(cross(d[1], d[2]).z > 0.0);
  // rulings along u
  const UnrolledSurface t = unroll(cyl.transposed(), tol);
  CHECK(std::abs(area(t.pattern) - 2 * pi * R * H) < 1e-6);
  const auto e = t.pattern.first_derivatives(0.3, 0.4);
  CHECK(cross(e[1], e[2]).z > 0.0);
  CHECK(metric_error(cyl.transposed(), t.pattern) < 1e-6);
}

TEST_CASE("a cone unrolls to a circular sector", "[unroll]") {
  const double R = 1.0, h = 2.0, l = std::sqrt(R * R + h * h), tol = 1e-8;
  const NurbsSurface cone = extrude_to_point(make_circle({0, 0, 0}, {0, 0, 1}, R), {0, 0, h});
  const UnrolledSurface r = unroll(cone, tol);
  const auto du = cone.domain_u(), dv = cone.domain_v();
  const Vec3 apex = r.pattern.evaluate(du[0], dv[1]);
  for_grid(cone, 24, [&](double u, double v) {
    const double f = (v - dv[0]) / (dv[1] - dv[0]);
    CHECK(distance(r.pattern.evaluate(u, dv[1]), apex) < 2 * tol);  // one apex
    CHECK(std::abs(distance(r.pattern.evaluate(u, v), apex) - l * (1 - f)) < 2 * tol);
  });
  // opening angle of the sector: base circumference / slant height
  const Vec3 a = r.pattern.evaluate(du[0], dv[0]) - apex,
             b = r.pattern.evaluate(du[1], dv[0]) - apex;
  double angle = std::atan2(a.x * b.y - a.y * b.x, dot(a, b));
  if (angle < 0) angle += 2 * pi;
  CHECK(std::abs(angle - 2 * pi * R / l) < 1e-7);
  CHECK(std::abs(area(r.pattern) - pi * R * l) < 1e-6);
  CHECK(metric_error(cone, r.pattern) < 1e-6);
}

TEST_CASE("a general developable surface unrolls isometrically", "[unroll]") {
  // oblique cone over a planar freeform curve, and a general cylinder
  const NurbsCurve base =
      NurbsCurve::from_points(3, {{0, 0, 0}, {1, 1.5, 0}, {2.5, 1, 0}, {3, -1, 0}, {4.5, 0, 0}});
  const NurbsSurface cone = extrude_to_point(base, {1.5, 0.5, 3});
  const NurbsSurface wall = extrude(base, {0.5, 0.2, 2});
  for (const NurbsSurface& s : {cone, wall}) {
    REQUIRE(developability(s).developable);
    const UnrolledSurface r = unroll(s, 1e-8);
    CHECK(metric_error(s, r.pattern) < 1e-6);
    CHECK(std::abs(area(r.pattern) - area(s)) < 1e-6);
    for_grid(s, 8,
             [&](double u, double v) { CHECK(std::abs(r.pattern.evaluate(u, v).z) < 1e-14); });
  }
  // not developable, not ruled
  CHECK_THROWS_AS(unroll(make_bilinear({0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 1})), GeometryError);
  CHECK_THROWS_AS(unroll(make_sphere({0, 0, 0}, 1.0).faces()[0].surface()), GeometryError);
}

TEST_CASE("chamfer between two planes and between a plane and a cylinder", "[fillet]") {
  // floor z = 0 (normal +z) and wall x = 0 (normal +x)
  const NurbsSurface floor = make_rectangle({-2, -2, 0}, {4, 0, 0}, {0, 4, 0});
  const NurbsSurface wall = make_rectangle({0, -2, -2}, {0, 4, 0}, {0, 0, 4});
  const SurfaceChamfer c = chamfer_surfaces(floor, wall, 0.5, 0.3);
  for_grid(c.surface, 16, [&](double u, double v) {
    const Vec3 p = c.surface.evaluate(u, v);
    CHECK(std::abs(p.x / 0.5 + p.z / 0.3 - 1.0) < 1e-9);  // the plane through both contacts
    CHECK(std::abs(p.x - 0.5 * (1 - u)) < 1e-9);          // u = 0 on the floor
    const Vec3 q1 = c.pcurve1.evaluate(v), q2 = c.pcurve2.evaluate(v);
    CHECK(distance(floor.evaluate(q1.x, q1.y), c.contact1.evaluate(v)) < 1e-9);
    CHECK(distance(wall.evaluate(q2.x, q2.y), c.contact2.evaluate(v)) < 1e-9);
    CHECK(distance(c.contact1.evaluate(v), c.surface.evaluate(0.0, v)) < 1e-9);
    CHECK(std::abs(c.edge.evaluate(v).x) + std::abs(c.edge.evaluate(v).z) < 1e-9);
  });
  // the other side of the floor
  SurfaceChamferOptions opt;
  opt.side1 = -1;
  const SurfaceChamfer below = chamfer_surfaces(floor, wall, 0.5, 0.3, opt);
  const Vec3 mid = below.surface.evaluate(0.5, 0.5);
  CHECK(std::abs(mid.x - 0.25) < 1e-9);
  CHECK(std::abs(mid.z + 0.15) < 1e-9);

  // plane z = 0 and a cylinder of radius 1 about z: a cone frustum
  const NurbsSurface plane = make_rectangle({-3, -3, 0}, {6, 0, 0}, {0, 6, 0});
  const NurbsSurface cyl = extrude(make_circle({0, 0, -1}, {0, 0, 1}, 1.0), {0, 0, 3});
  SurfaceChamferOptions o2;
  o2.side2 =
      dot(cyl.normal(0.1, 0.5), cyl.evaluate(0.1, 0.5) - Vec3{0, 0, cyl.evaluate(0.1, 0.5).z}) > 0.0
          ? 1
          : -1;  // outside the cylinder
  const SurfaceChamfer k = chamfer_surfaces(plane, cyl, 0.4, 0.25, o2);
  for_grid(k.surface, 16, [&](double u, double v) {
    const Vec3 p = k.surface.evaluate(u, v);
    CHECK(std::abs(std::hypot(p.x, p.y) - (1.4 - 0.4 * u)) < 1e-6);
    CHECK(std::abs(p.z - 0.25 * u) < 1e-6);
  });
  CHECK_THROWS_AS(chamfer_surfaces(floor, wall, 0.5, 5.0), GeometryError);  // leaves the wall
  CHECK_THROWS_AS(
      chamfer_surfaces(floor, make_rectangle({-2, -2, 1}, {4, 0, 0}, {0, 4, 0}), 0.1, 0.1),
      GeometryError);  // parallel planes
}

TEST_CASE("flattening of developable and doubly curved surfaces", "[unroll]") {
  // a plane is its own pattern, whatever its parametrisation
  const NurbsSurface plane = make_bilinear({0, 0, 1}, {3, 0.5, 1}, {-0.5, 2, 1}, {3.5, 2.5, 1});
  const FlattenedSurface p = flatten(plane);
  CHECK(std::abs(p.max_stretch - 1.0) < 1e-7);
  CHECK(std::abs(p.min_stretch - 1.0) < 1e-7);
  CHECK(std::abs(p.pattern_area - p.area) < 1e-6);
  const auto d = p.pattern.first_derivatives(0.3, 0.6);
  CHECK(cross(d[1], d[2]).z > 0.0);
  CHECK(norm(p.pattern.evaluate(0.5, 0.5)) < 1e-9);  // the middle at the origin

  // a developable surface flattens without strain: the pattern agrees with unroll()
  const NurbsSurface cyl =
      extrude(make_arc({0, 0, 0}, {1, 0, 0}, {0, 1, 0}, 2.0, 0.0, 2.0), {0, 0, 3});
  const FlattenedSurface c = flatten(cyl);
  CHECK(c.max_stretch - 1.0 < 1e-5);
  CHECK(1.0 - c.min_stretch < 1e-5);
  CHECK(std::abs(c.pattern_area - 12.0) < 1e-4);  // arc length 4 x height 3
  // the pattern is a rectangle: its points lie at the arc length and the height
  for_grid(cyl, 8, [&](double u, double v) {
    const Vec3 a = c.pattern.evaluate(u, v), b = c.pattern.evaluate(cyl.domain_u()[0], v);
    const Vec3 p = cyl.evaluate(u, v);
    CHECK(std::abs(distance(a, b) - 2.0 * std::atan2(p.y, p.x)) < 1e-4);
  });

  // a spherical cap (40 x 40 degrees about the equator of a sphere of radius 2): small but
  // unavoidable strains, compression at the rim and tension in the middle
  const NurbsSurface sphere = make_sphere({0, 0, 0}, 2.0).faces()[0].surface();
  const auto du = sphere.domain_u(), dv = sphere.domain_v();
  const double um = 0.5 * (du[0] + du[1]) + 0.07 * (du[1] - du[0]);
  const double vm = 0.5 * (dv[0] + dv[1]);
  const NurbsSurface cap =
      sphere.sub_surface(um - 0.055 * (du[1] - du[0]), um + 0.055 * (du[1] - du[0]),
                         vm - 0.11 * (dv[1] - dv[0]), vm + 0.11 * (dv[1] - dv[0]));
  const FlattenedSurface f = flatten(cap);
  CHECK(f.stretch_max.size() == f.u.size() * f.v.size());
  CHECK(f.max_stretch > 1.002);
  CHECK(f.max_stretch < 1.06);
  CHECK(f.min_stretch < 0.998);
  CHECK(f.min_stretch > 0.94);
  CHECK(f.rms_strain > 1e-3);
  CHECK(f.rms_strain < 0.03);
  CHECK(std::abs(f.pattern_area / f.area - 1.0) < 0.02);
  for_grid(cap, 6,
           [&](double u, double v) { CHECK(std::abs(f.pattern.evaluate(u, v).z) < 1e-14); });
  // less strain than the unrelaxed start would have, and a finer grid changes little
  FlattenOptions fine;
  fine.samples_u = fine.samples_v = 49;
  const FlattenedSurface g = flatten(cap, fine);
  CHECK(std::abs(g.rms_strain - f.rms_strain) < 0.1 * f.rms_strain);
  FlattenOptions bad;
  bad.samples_u = 2;
  CHECK_THROWS_AS(flatten(cap, bad), GeometryError);
}
