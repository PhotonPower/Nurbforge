#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <numbers>
#include <vector>

#include "nurbforge/analysis/closest_point.hpp"
#include "nurbforge/core/error.hpp"
#include "nurbforge/fitting/fit.hpp"
#include "nurbforge/geometry/primitives.hpp"
#include "nurbforge/modeling/blend.hpp"
#include "nurbforge/modeling/construct.hpp"

using namespace nf;
constexpr double pi = std::numbers::pi;

namespace {

// Bicubic surface through z = f(x, y) on [x0, x0 + w] x [y0, y0 + h].
NurbsSurface wavy(double x0, double y0, double w, double h, double z0, double phase) {
  std::vector<Vec3> pts;
  const std::size_t n = 7;
  for (std::size_t i = 0; i < n; ++i)
    for (std::size_t j = 0; j < n; ++j) {
      const double x = x0 + w * static_cast<double>(i) / (n - 1);
      const double y = y0 + h * static_cast<double>(j) / (n - 1);
      pts.push_back({x, y, z0 + 0.3 * std::sin(1.3 * x + phase) * std::cos(0.9 * y)});
    }
  return interpolate_surface(pts, n, n, 3, 3);
}

// Normal curvature of s at (u, v) in the tangent direction t (normal oriented along n_ref).
double normal_curvature(const NurbsSurface& s, double u, double v, const Vec3& t,
                        const Vec3& n_ref) {
  const auto d = s.derivatives(u, v, 2);
  Vec3 n = normalized(cross(d[1][0], d[0][1]));
  if (dot(n, n_ref) < 0) n = -1.0 * n;
  const double E = dot(d[1][0], d[1][0]), F = dot(d[1][0], d[0][1]), G = dot(d[0][1], d[0][1]);
  const double ta = dot(t, d[1][0]), tb = dot(t, d[0][1]);
  const double det = E * G - F * F;
  const double a = (G * ta - F * tb) / det, b = (E * tb - F * ta) / det;  // t = a S_u + b S_v
  const double L = dot(d[2][0], n), M = dot(d[1][1], n), N = dot(d[0][2], n);
  return (L * a * a + 2 * M * a * b + N * b * b) / (E * a * a + 2 * F * a * b + G * b * b);
}

NurbsSurface reversed_v(const NurbsSurface& s) {
  const auto& kv = s.knots_v();
  KnotVector k(kv.size());
  for (std::size_t i = 0; i < kv.size(); ++i) k[i] = kv.front() + kv.back() - kv[kv.size() - 1 - i];
  std::vector<Vec3> P;
  std::vector<double> W;
  for (std::size_t i = 0; i < s.num_u(); ++i)
    for (std::size_t j = s.num_v(); j-- > 0;) {
      P.push_back(s.control_point(i, j));
      W.push_back(s.weight(i, j));
    }
  return {s.degree_u(), s.degree_v(), s.knots_u(), k, s.num_u(), s.num_v(), P, W};
}

struct JoinError {
  double gap{0}, angle{0}, curvature{0};
};

// Deviations between surface b along its iso-line v = vb (b's edge) and the surface s.
JoinError join_error(const NurbsSurface& b, double vb, const NurbsSurface& s) {
  JoinError e;
  const auto [u0, u1] = b.domain_u();
  for (int k = 1; k < 20; ++k) {
    const double u = u0 + (u1 - u0) * k / 20.0;
    const auto d = b.derivatives(u, vb, 1);
    const SurfacePoint f = closest_point(s, d[0][0]);
    e.gap = std::max(e.gap, f.distance);
    const Vec3 nb = normalized(cross(d[1][0], d[0][1]));
    const Vec3 ns = s.normal(f.u, f.v);
    e.angle = std::max(e.angle, norm(cross(nb, ns)));
    const Vec3 along = normalized(d[1][0]);
    const Vec3 across = normalized(cross(nb, along));
    for (const Vec3& t : {across, normalized(across + along), normalized(across - 0.5 * along)}) {
      const double kb = normal_curvature(b, u, vb, t, nb);
      const double ks = normal_curvature(s, f.u, f.v, t, nb);
      e.curvature = std::max(e.curvature, std::abs(kb - ks));
    }
  }
  return e;
}

}  // namespace

TEST_CASE("blend surfaces: G0, G1, G2 transitions", "[blend][surface]") {
  const NurbsSurface a = wavy(0, 0, 2, 3, 0, 0.0);
  const NurbsSurface b = wavy(3, -0.5, 2, 4, 0.8, 1.0);

  SECTION("G2 between two bicubic surfaces (opposite orientation handled)") {
    for (const bool flip : {false, true}) {
      const NurbsSurface bb = flip ? reversed_v(b) : b;
      const auto bl = blend_surfaces(a, SurfaceSide::UMax, bb, SurfaceSide::UMin, 2);
      CHECK(bl.degree_v() == 5);
      const auto e0 = join_error(bl, 0.0, a), e1 = join_error(bl, 1.0, bb);
      CHECK(e0.gap <= 1e-9);
      CHECK(e1.gap <= 1e-9);
      CHECK(e0.angle <= 1e-9);
      CHECK(e1.angle <= 1e-9);
      CHECK(e0.curvature <= 1e-7);
      CHECK(e1.curvature <= 1e-7);
      CHECK(distance(
                bl.evaluate(0.5, 0.5),
                blend_surfaces(a, SurfaceSide::UMax, b, SurfaceSide::UMin, 2).evaluate(0.5, 0.5)) <=
            1e-12);
    }
  }
  SECTION("G1: cubic, tangent planes agree, curvature does not") {
    const auto bl = blend_surfaces(a, SurfaceSide::UMax, b, SurfaceSide::UMin, 1);
    CHECK(bl.degree_v() == 3);
    const auto e0 = join_error(bl, 0.0, a), e1 = join_error(bl, 1.0, b);
    CHECK(e0.angle <= 1e-9);
    CHECK(e1.angle <= 1e-9);
    CHECK(e0.curvature > 1e-3);
  }
  SECTION("G0: ruled") {
    const auto bl = blend_surfaces(a, SurfaceSide::UMax, b, SurfaceSide::UMin, 0);
    CHECK(bl.degree_v() == 1);
    CHECK(join_error(bl, 0.0, a).gap <= 1e-9);
    CHECK(join_error(bl, 1.0, b).gap <= 1e-9);
  }
  SECTION("exact for rational sides with constant weights across: cylinder rim to a plane") {
    // same circle parametrisation on both sides -> a rotationally symmetric blend
    const NurbsSurface cyl = revolve(make_line({1, 0, 0}, {1, 0, 2}), {0, 0, 0}, {0, 0, 1}, 2 * pi);
    const NurbsSurface ring = revolve(make_line({2, 0, 3}, {3, 0, 3}), {0, 0, 0}, {0, 0, 1},
                                      2 * pi);  // annulus at z = 3, v radial
    const auto bl = blend_surfaces(cyl, SurfaceSide::VMax, ring, SurfaceSide::VMin, 2);
    CHECK(!bl.weights().empty());
    const auto e0 = join_error(bl, 0.0, cyl), e1 = join_error(bl, 1.0, ring);
    CHECK(e0.gap <= 1e-9);
    CHECK(e1.gap <= 1e-9);
    CHECK(e0.angle <= 1e-9);
    CHECK(e1.angle <= 1e-9);
    CHECK(e0.curvature <= 1e-7);
    CHECK(e1.curvature <= 1e-7);
    // rotational symmetry: every cross-section is the same curve rotated
    const Vec3 p = bl.evaluate(0.0, 0.5), q = bl.evaluate(0.37, 0.5);
    CHECK(std::abs(std::hypot(p.x, p.y) - std::hypot(q.x, q.y)) <= 1e-12);
    CHECK(std::abs(p.z - q.z) <= 1e-12);
  }
  SECTION("sampled rows when the weights vary across: a sphere zone") {
    const NurbsSurface zone = revolve(make_arc({0, 0, 0}, {1, 0, 0}, {0, 0, 1}, 1.0, 0, pi / 3),
                                      {0, 0, 0}, {0, 0, 1}, 2 * pi);
    const NurbsSurface top =
        revolve(make_line({0.2, 0, 1.6}, {0.05, 0, 1.6}), {0, 0, 0}, {0, 0, 1}, 2 * pi);
    const auto bl = blend_surfaces(zone, SurfaceSide::VMax, top, SurfaceSide::VMin, 2, 1.0, 1e-8);
    const auto e0 = join_error(bl, 0.0, zone), e1 = join_error(bl, 1.0, top);
    CHECK(e0.gap <= 1e-7);
    CHECK(e1.gap <= 1e-9);
    CHECK(e0.angle <= 1e-6);
    CHECK(e0.curvature <= 2e-4);
    CHECK(e1.curvature <= 1e-7);
  }
  SECTION("errors") {
    CHECK_THROWS_AS(blend_surfaces(a, SurfaceSide::UMax, b, SurfaceSide::UMin, 3), GeometryError);
    CHECK_THROWS_AS(blend_surfaces(a, SurfaceSide::UMax, a, SurfaceSide::UMax), GeometryError);
    CHECK_THROWS_AS(blend_surfaces(a, SurfaceSide::UMax, b, SurfaceSide::UMin, 2, 0.0),
                    GeometryError);
  }
}

TEST_CASE("match surface: G0, G1, G2 to a target side", "[blend][surface]") {
  const NurbsSurface target = wavy(0, 0, 2, 3, 0, 0.0);
  // a surface whose VMin side lies near the target's UMax side (x = 2)
  const NurbsSurface s = wavy(2.3, -0.2, 2, 3.4, 0.4, 2.0).transposed();

  for (int c = 0; c <= 2; ++c) {
    const auto m = match_surface(s, SurfaceSide::VMin, target, SurfaceSide::UMax, c);
    CHECK(m.domain_u() == s.domain_u());
    CHECK(m.domain_v() == s.domain_v());
    const auto e = join_error(m, m.domain_v()[0], target);
    CHECK(e.gap <= 1e-9);
    if (c >= 1) CHECK(e.angle <= 1e-9);
    if (c == 2) CHECK(e.curvature <= 1e-7);
    if (c < 2) CHECK(e.curvature > 1e-4);
    // the opposite side and its cross derivatives up to order c are kept
    double dev = 0.0;
    const double v1 = s.domain_v()[1];
    for (int k = 0; k <= 10; ++k) {
      const double u = k / 10.0;
      const auto d0 = s.derivatives(u, v1, 2), d1 = m.derivatives(u, v1, 2);
      for (int o = 0; o <= c; ++o) dev = std::max(dev, distance(d0[0][o], d1[0][o]));
    }
    CHECK(dev <= 1e-9);
  }
  SECTION("all four sides of the surface") {
    const NurbsSurface sq = wavy(2.1, 0, 2, 3, 0.3, 0.5);
    for (const auto side :
         {SurfaceSide::UMin, SurfaceSide::UMax, SurfaceSide::VMin, SurfaceSide::VMax}) {
      const auto m = match_surface(sq, side, target, SurfaceSide::UMax, 1);
      const double v = side == SurfaceSide::VMin ? m.domain_v()[0] : m.domain_v()[1];
      const NurbsSurface mm =
          side == SurfaceSide::UMin || side == SurfaceSide::UMax ? m.transposed() : m;
      const double vv = side == SurfaceSide::UMin   ? mm.domain_v()[0]
                        : side == SurfaceSide::UMax ? mm.domain_v()[1]
                                                    : v;
      const auto e = join_error(mm, vv, target);
      CHECK(e.gap <= 1e-9);
      CHECK(e.angle <= 1e-9);
    }
  }
}
