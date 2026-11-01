#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <numbers>
#include <vector>

#include "nurbforge/core/error.hpp"
#include "nurbforge/geometry/primitives.hpp"
#include "nurbforge/modeling/sweep.hpp"

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

NurbsCurve space_curve() {
  return NurbsCurve::from_points(3, {{0, 0, 0}, {2, 0, 1}, {3, 2, 2}, {2, 4, 2}, {0, 5, 3}, {-1, 4, 5}});
}

}  // namespace

TEST_CASE("translational sweep is exact", "[sweep]") {
  const auto profile = make_circle({0, 0, 0}, {0, 0, 1}, 0.5);  // rational
  const auto path = space_curve();
  const auto s = sweep_translational(profile, path);
  const Vec3 start = path.evaluate(0.0);
  for_grid(s, 20, [&](double u, double v) {
    CHECK(distance(s.evaluate(u, v), profile.evaluate(u) + path.evaluate(v) - start) < 1e-12);
  });
  // SweepFrame::Fixed without twist returns the exact translational sweep
  const auto f = sweep(profile, path, {SweepFrame::Fixed});
  CHECK(distance(f.evaluate(0.3, 0.7), s.evaluate(0.3, 0.7)) < 1e-14);
}

TEST_CASE("circle swept along a straight line is a cylinder", "[sweep]") {
  const auto s = sweep(make_circle({0, 0, 0}, {0, 0, 1}, 2.0), make_line({0, 0, 0}, {0, 0, 5}));
  for_grid(s, 20, [&](double u, double v) {
    const Vec3 p = s.evaluate(u, v);
    CHECK(std::abs(std::hypot(p.x, p.y) - 2.0) < 1e-12);
    CHECK(std::abs(p.z - 5.0 * v) < 1e-12);
  });
}

TEST_CASE("circle swept along a circular arc is a torus within the tolerance", "[sweep]") {
  const double R = 5.0, r = 1.0;
  const auto path = make_arc({0, 0, 0}, {1, 0, 0}, {0, 1, 0}, R, 0.0, 3 * pi / 2);
  const auto profile = make_circle({R, 0, 0}, {0, 1, 0}, r);  // in the plane normal to the path
  for (double tol : {1e-5, 1e-8}) {
    SweepOptions opt;
    opt.tolerance = tol;
    const auto s = sweep(profile, path, opt);
    double worst = 0.0;
    for_grid(s, 60, [&](double u, double v) {
      const Vec3 p = s.evaluate(u, v);
      worst = std::max(worst, std::abs(std::hypot(std::hypot(p.x, p.y) - R, p.z) - r));
    });
    CHECK(worst < 5 * tol);
  }
}

TEST_CASE("rotation-minimising frame keeps a planar wall upright", "[sweep]") {
  const auto path = NurbsCurve::from_points(3, {{0, 0, 0}, {3, 0, 0}, {3, 3, 0}, {0, 4, 0}, {-2, 1, 0}});
  const auto profile = make_line({0, 0, 0}, {0, 0, 1});  // normal to the plane of the path
  const auto s = sweep(profile, path);
  for_grid(s, 30, [&](double u, double v) {
    const Vec3 p = s.evaluate(u, v), q = path.evaluate(v);
    CHECK(std::abs(p.z - u) < 1e-9);
    CHECK(std::hypot(p.x - q.x, p.y - q.y) < 1e-9);
  });
}

TEST_CASE("rotation-minimising frame has no spurious twist on a space curve", "[sweep]") {
  const auto path = space_curve();
  const Vec3 t0 = normalized(path.derivatives(0.0, 1)[1]);
  const Vec3 n0 = normalized(cross(t0, {0, 0, 1}));
  SweepOptions opt;
  opt.tolerance = 1e-9;
  const auto s = sweep(make_line({0, 0, 0}, n0), path, opt);  // S(u, v) = T(v) + u r(v)
  for (double v = 0.05; v < 1.0; v += 0.1) {
    const auto D = s.derivatives(0.0, v, 2);
    const Vec3 r = D[1][0];            // transported vector, |r| = 1
    const Vec3 dr = D[1][1];           // its derivative along the path
    const Vec3 t = normalized(path.derivatives(v, 1)[1]);
    CHECK(std::abs(norm(r) - 1.0) < 1e-6);
    CHECK(std::abs(dot(r, t)) < 1e-6);
    CHECK(std::abs(dot(dr, cross(t, r))) < 1e-4 * (1.0 + norm(dr)));  // r' has no binormal part
  }
}

TEST_CASE("sweep with twist", "[sweep]") {
  SweepOptions opt;
  opt.twist = pi / 2;
  const auto s = sweep(make_line({0, 0, 0}, {1, 0, 0}), make_line({0, 0, 0}, {0, 0, 5}), opt);
  CHECK(distance(s.evaluate(1.0, 1.0), {0, 1, 5}) < 1e-6);
  CHECK(distance(s.evaluate(1.0, 0.5), {std::cos(pi / 4), std::sin(pi / 4), 2.5}) < 1e-6);
}

TEST_CASE("two-rail sweep scales the profile between the rails", "[sweep]") {
  const auto rail1 = make_line({0, 0, 0}, {0, 10, 0});
  const auto rail2 = make_line({1, 0, 0}, {3, 10, 0});
  // half circle from rail1(0) to rail2(0), bulging upwards (+z)
  const auto profile = make_arc({0.5, 0, 0}, {-1, 0, 0}, {0, 0, 1}, 0.5, 0.0, pi);
  const auto s = sweep2(profile, rail1, rail2);
  for (double v = 0.0; v <= 1.0; v += 0.1) {
    const Vec3 a = rail1.evaluate(v), b = rail2.evaluate(v);
    CHECK(distance(s.evaluate(0.0, v), a) < 1e-6);
    CHECK(distance(s.evaluate(1.0, v), b) < 1e-6);
    // the top of the arc lies above the midpoint at half the rail distance
    const Vec3 top = s.evaluate(0.5, v);
    CHECK(distance(top, 0.5 * (a + b) + Vec3{0, 0, 0.5 * distance(a, b)}) < 1e-6);
  }
}

TEST_CASE("sweep reports an unreachable tolerance", "[sweep]") {
  SweepOptions opt;
  opt.tolerance = 1e-14;
  opt.max_sections = 9;
  const auto path = make_arc({0, 0, 0}, {1, 0, 0}, {0, 1, 0}, 5.0, 0.0, pi);
  CHECK_THROWS_AS(sweep(make_circle({5, 0, 0}, {0, 1, 0}, 1.0), path, opt), GeometryError);
}
