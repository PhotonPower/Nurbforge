#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <numbers>
#include <thread>
#include <vector>

#include "nurbforge/analysis/shell_section.hpp"
#include "nurbforge/core/error.hpp"
#include "nurbforge/geometry/primitives.hpp"
#include "nurbforge/modeling/construct.hpp"
#include "nurbforge/topology/boolean.hpp"

using namespace nf;
using Catch::Approx;
constexpr double pi = std::numbers::pi;

namespace {

Shell box_shell(double h) {  // cube [-h, h]^3
  const double s = 2 * h;
  return Shell::sew({TrimmedSurface(make_rectangle({-h, -h, -h}, {s, 0, 0}, {0, s, 0})),
                     TrimmedSurface(make_rectangle({-h, -h, -h}, {s, 0, 0}, {0, 0, s})),
                     TrimmedSurface(make_rectangle({-h, h, -h}, {s, 0, 0}, {0, 0, s})),
                     TrimmedSurface(make_rectangle({-h, -h, -h}, {0, s, 0}, {0, 0, s})),
                     TrimmedSurface(make_rectangle({h, -h, -h}, {0, s, 0}, {0, 0, s})),
                     TrimmedSurface(make_rectangle({-h, -h, h}, {s, 0, 0}, {0, s, 0}))});
}

Shell sphere_shell(const Vec3& c, double R) {
  return Shell::sew({TrimmedSurface(
      revolve(make_arc(c, {0, 0, -1}, {1, 0, 0}, R, 0.0, pi), c, {0, 0, 1}, 2 * pi))});
}

Shell cylinder_shell(double R, double z0, double z1) {  // axis z
  const auto side = extrude(make_circle({0, 0, z0}, {0, 0, 1}, R), {0, 0, z1 - z0});
  const auto lid = [&](double z) {
    return TrimmedSurface(make_rectangle({-R, -R, z}, {2 * R, 0, 0}, {0, 2 * R, 0}),
                          circle_loop(0.5, 0.5, 0.5));
  };
  return Shell::sew({TrimmedSurface(side), lid(z0), lid(z1)});
}

Shell frustum_shell(double r0, double r1, double h) {  // radius r0 at z = 0, r1 at z = h
  const auto side = revolve(make_line({r0, 0, 0}, {r1, 0, h}), {0, 0, 0}, {0, 0, 1}, 2 * pi);
  const auto lid = [&](double z, double r) {
    return TrimmedSurface(make_rectangle({-r, -r, z}, {2 * r, 0, 0}, {0, 2 * r, 0}),
                          circle_loop(0.5, 0.5, 0.5));
  };
  return Shell::sew({TrimmedSurface(side), lid(0.0, r0), lid(h, r1)});
}

Shell torus_shell(double R, double r, const Vec3& axis) {  // centre 0
  const Vec3 a = normalized(axis);
  const Vec3 x =
      std::abs(a.x) < 0.9 ? normalized(cross(a, {1, 0, 0})) : normalized(cross(a, {0, 1, 0}));
  const Vec3 y = cross(a, x);
  // profile circle in the (x, a) plane about R x, revolved about a
  return Shell::sew({TrimmedSurface(revolve(make_circle(R * x, y, r), {0, 0, 0}, a, 2 * pi))});
}

double total_area(const PlanarSection& s) {
  double a = 0;
  for (const auto& l : s.loops) a += l.area;
  return a;
}

}  // namespace

TEST_CASE("shell section: box and cylinder (exact pieces)", "[shell_section]") {
  const auto box = section(box_shell(1.0), {{0, 0, 0.3}, {0, 0, 1}});
  REQUIRE(box.loops.size() == 1);
  CHECK(box.loops[0].area == Approx(4.0).epsilon(1e-14));
  CHECK(box.loops[0].curves.size() == 4);
  for (const auto& c : box.loops[0].curves) CHECK(c.exact);
  const auto cyl = section(cylinder_shell(1.0, 0.0, 2.0), {{0, 0, 0.7}, {0, 0, 1}});
  REQUIRE(cyl.loops.size() == 1);
  CHECK(cyl.loops[0].area == Approx(pi).epsilon(1e-13));
  CHECK(cyl.loops[0].curves.front().exact);
  // every point of the section lies on the circle, the normal is radial and horizontal
  for (const auto& c : cyl.loops[0].curves)
    for (double t : {0.1, 0.5, 0.9}) {
      const double u = c.curve.domain_start() + t * (c.curve.domain_end() - c.curve.domain_start());
      const Vec3 p = c.curve.evaluate(u), q = c.curve3d.evaluate(u), n = c.normal(u);
      CHECK(std::hypot(p.x, p.y) == Approx(1.0).epsilon(1e-14));
      CHECK(q.z == Approx(0.7).epsilon(1e-14));
      CHECK(dot(n, normalized(Vec3{q.x, q.y, 0})) == Approx(1.0).epsilon(1e-12));
    }
}

TEST_CASE("shell section: sphere, frustum, torus against analytic areas", "[shell_section]") {
  const auto sph = section(sphere_shell({0, 0, 0}, 1.0), {{0, 0, 0.5}, {0, 0, 1}});
  REQUIRE(sph.loops.size() == 1);
  CHECK(sph.loops[0].area == Approx(0.75 * pi).epsilon(1e-12));
  const auto fr = section(frustum_shell(1.0, 0.4, 2.0), {{0, 0, 0.5}, {0, 0, 1}});
  REQUIRE(fr.loops.size() == 1);
  CHECK(fr.loops[0].area == Approx(pi * 0.85 * 0.85).epsilon(1e-12));
  // frustum wall: tilted outward normal (cos t, sin t, (r0 - r1)/h) / |.|
  const auto& c = fr.loops[0].curves.front();
  const double u = 0.5 * (c.curve.domain_start() + c.curve.domain_end());
  CHECK(c.normal(u).z == Approx(0.3 / std::sqrt(1.09)).epsilon(1e-12));
  // torus R = 2, r = 0.5 about z cut at z = 0.3: outer and inner loop
  const auto tor = section(torus_shell(2.0, 0.5, {0, 0, 1}), {{0, 0, 0.3}, {0, 0, 1}});
  REQUIRE(tor.loops.size() == 2);
  const double rho = std::sqrt(0.25 - 0.09);
  const double outer = std::max(tor.loops[0].area, tor.loops[1].area),
               inner = std::min(tor.loops[0].area, tor.loops[1].area);
  CHECK(outer == Approx(pi * (2 + rho) * (2 + rho)).epsilon(1e-12));
  CHECK(inner == Approx(-pi * (2 - rho) * (2 - rho)).epsilon(1e-12));
}

TEST_CASE("shell section: box with a spherical cavity (hole), tilted plane (traced)",
          "[shell_section]") {
  const Shell cav = boolean(box_shell(1.0), sphere_shell({0, 0, 0}, 0.5), BooleanOp::Difference);
  const auto s = section(cav, {{0, 0, 0.2}, {0, 0, 1}});
  REQUIRE(s.loops.size() == 2);
  CHECK(total_area(s) == Approx(4.0 - pi * (0.25 - 0.04)).epsilon(1e-9));
  CHECK(std::min(s.loops[0].area, s.loops[1].area) == Approx(-pi * 0.21).epsilon(1e-9));
  // a plane at distance 0.3 with a tilted normal: the sphere face is traced (not isoparametric)
  const Vec3 n = normalized(Vec3{1, 0, 1});
  const auto t = section(sphere_shell({0, 0, 0}, 1.0), {0.3 * n, n});
  REQUIRE(t.loops.size() == 1);
  CHECK_FALSE(t.loops[0].curves.front().exact);
  CHECK(t.loops[0].area == Approx(pi * (1 - 0.09)).epsilon(1e-8));
  // a plane containing a face is rejected
  CHECK_THROWS_AS(section(box_shell(1.0), {{0, 0, 1.0}, {0, 0, 1}}), GeometryError);
}

TEST_CASE("critical heights and extrusion bands", "[shell_section]") {
  const auto hs = critical_heights(sphere_shell({0, 0, 0}, 1.0), {0, 0, 1}, {1e-9});
  REQUIRE(hs.size() == 2);
  CHECK(hs[0] == Approx(-1.0).margin(1e-12));
  CHECK(hs[1] == Approx(1.0).margin(1e-12));
  // torus with its axis perpendicular to the direction: four heights -(R + r), -(R - r), R - r, R +
  // r
  const auto ht = critical_heights(torus_shell(2.0, 0.5, {1, 0, 0}), {0, 0, 1}, {1e-9});
  REQUIRE(ht.size() == 4);
  CHECK(ht[0] == Approx(-2.5).margin(1e-9));
  CHECK(ht[1] == Approx(-1.5).margin(1e-9));
  CHECK(ht[2] == Approx(1.5).margin(1e-9));
  CHECK(ht[3] == Approx(2.5).margin(1e-9));
  // torus about z: the top and bottom circles
  const auto hz = critical_heights(torus_shell(2.0, 0.5, {0, 0, 1}), {0, 0, 1}, {1e-9});
  REQUIRE(hz.size() == 2);
  CHECK(hz[0] == Approx(-0.5).margin(1e-12));
  // cylinder: the lids
  const auto hc = critical_heights(cylinder_shell(1.0, 0.0, 2.0), {0, 0, 1}, {1e-9});
  REQUIRE(hc.size() == 2);
  CHECK(hc[1] == Approx(2.0).margin(1e-12));
  CHECK(is_extrusion_band(cylinder_shell(1.0, 0.0, 2.0), {0, 0, 1}, 0.0, 2.0));
  CHECK(is_extrusion_band(box_shell(1.0), {0, 0, 1}, -1.0, 1.0));
  CHECK_FALSE(is_extrusion_band(frustum_shell(1.0, 0.4, 2.0), {0, 0, 1}, 0.0, 2.0));
  CHECK_FALSE(is_extrusion_band(sphere_shell({0, 0, 0}, 1.0), {0, 0, 1}, -0.5, 0.5));
  // a band outside the body: trivially an extrusion (nothing to check)
  CHECK(is_extrusion_band(sphere_shell({0, 0, 0}, 1.0), {0, 0, 1}, 2.0, 3.0));
}

TEST_CASE("shell section: concurrent const queries", "[shell_section][threads]") {
  // all section queries are const and thread-safe: same results from four threads
  const Shell cav = boolean(box_shell(1.0), sphere_shell({0, 0, 0}, 0.5), BooleanOp::Difference);
  const double ref = total_area(section(cav, {{0, 0, 0.1}, {0, 0, 1}}));
  std::vector<double> got(4, 0.0);
  std::vector<std::thread> th;
  for (std::size_t k = 0; k < got.size(); ++k)
    th.emplace_back([&, k] {
      got[k] = total_area(section(cav, {{0, 0, 0.1}, {0, 0, 1}}));
      (void)critical_heights(cav, {0, 0, 1});
    });
  for (auto& t : th) t.join();
  for (double g : got) CHECK(g == ref);
}
