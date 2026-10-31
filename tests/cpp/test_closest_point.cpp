#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numbers>
#include <vector>

#include "nurbforge/analysis/closest_point.hpp"
#include "nurbforge/geometry/primitives.hpp"
#include "nurbforge/modeling/construct.hpp"

using namespace nf;

namespace {

// Deterministic pseudo-random numbers in [lo, hi].
struct Rng {
  std::uint64_t s = 12345;
  double operator()(double lo, double hi) {
    s = s * 6364136223846793005ULL + 1442695040888963407ULL;
    return lo + (hi - lo) * static_cast<double>(s >> 11) / static_cast<double>(1ULL << 53);
  }
};

NurbsCurve rational_cubic() {
  return {3,
          {0, 0, 0, 0, 0.3, 0.3, 0.7, 1, 1, 1, 1},
          {{0, 0, 0}, {1, 2, 0}, {2, 3, 1}, {3, 1, 1}, {4, 0, 2}, {5, 2, 2}, {6, 1, 0}},
          {1.0, 0.5, 2.0, 1.0, 0.7, 1.3, 1.0}};
}

NurbsCurve wavy_cubic() {
  return NurbsCurve::from_points(
      3, {{0, 0, 0}, {1, 3, 0}, {2, -2, 0}, {3, 3, 1}, {4, -1, 0}, {5, 2, -1}, {6, 0, 0}});
}

NurbsSurface rational_patch() {
  const std::size_t nu = 5, nv = 4;
  std::vector<Vec3> P;
  std::vector<double> W;
  for (std::size_t i = 0; i < nu; ++i)
    for (std::size_t j = 0; j < nv; ++j) {
      const double x = static_cast<double>(i), y = static_cast<double>(j);
      P.push_back({x, y, 1.5 * std::sin(1.4 * x) * std::cos(1.3 * y)});
      W.push_back(1.0 + 0.4 * std::sin(1.7 * x + 0.9 * y) * std::sin(1.7 * x + 0.9 * y));
    }
  return {3, 2, {0, 0, 0, 0, 0.4, 1, 1, 1, 1}, {0, 0, 0, 0.6, 1, 1, 1}, nu, nv, P, W};
}

double brute_force(const NurbsCurve& c, const Vec3& P, std::size_t n = 200000) {
  double d = std::numeric_limits<double>::max();
  const double a = c.domain_start(), b = c.domain_end();
  for (std::size_t i = 0; i <= n; ++i)
    d = std::min(d, distance(c.evaluate(a + (b - a) * static_cast<double>(i) / static_cast<double>(n)), P));
  return d;
}

double brute_force(const NurbsSurface& s, const Vec3& P, std::size_t n = 400) {
  double d = std::numeric_limits<double>::max();
  const auto du = s.domain_u(), dv = s.domain_v();
  for (std::size_t i = 0; i <= n; ++i)
    for (std::size_t j = 0; j <= n; ++j) {
      const double u = du[0] + (du[1] - du[0]) * static_cast<double>(i) / static_cast<double>(n);
      const double v = dv[0] + (dv[1] - dv[0]) * static_cast<double>(j) / static_cast<double>(n);
      d = std::min(d, distance(s.evaluate(u, v), P));
    }
  return d;
}

}  // namespace

// ------------------------------------------------------------------ curves

TEST_CASE("curve point inversion", "[closest]") {
  for (const auto& c : {rational_cubic(), wavy_cubic()}) {
    for (double u = 0.0; u <= 1.0; u += 0.037) {
      const auto r = closest_point(c, c.evaluate(u));
      CHECK(r.distance < 1e-10);
      CHECK(std::abs(r.u - u) < 1e-8);
      CHECK(distance(r.point, c.evaluate(r.u)) < 1e-14);
    }
  }
}

TEST_CASE("curve closest point is never worse than brute force", "[closest]") {
  Rng rng;
  const NurbsCurve unclamped{2, {0, 1, 2, 3, 4, 5, 6, 7}, {{0, 0, 0}, {1, 1, 0}, {2, -1, 0}, {3, 1, 1}, {4, 0, 0}}, {}};
  for (const auto& c : {rational_cubic(), wavy_cubic(), unclamped}) {
    for (int k = 0; k < 25; ++k) {
      const Vec3 P{rng(-1, 7), rng(-3, 4), rng(-2, 3)};
      const auto r = closest_point(c, P);
      CHECK(r.distance <= brute_force(c, P) + 1e-9);
      CHECK(std::abs(distance(c.evaluate(r.u), P) - r.distance) < 1e-14);
      // first-order optimality in the interior: C' orthogonal to C - P
      if (r.u > c.domain_start() + 1e-9 && r.u < c.domain_end() - 1e-9 && r.distance > 1e-9) {
        const Vec3 t = c.derivatives(r.u, 1)[1];
        CHECK(std::abs(dot(t, r.point - P)) <= 1e-7 * norm(t) * r.distance);
      }
    }
  }
}

TEST_CASE("closest point on an exact circle", "[closest]") {
  const Vec3 C{1, -2, 0.5};
  const double R = 3.0;
  const auto circle = make_circle(C, {0, 0, 1}, R);
  Rng rng;
  for (int k = 0; k < 40; ++k) {
    const Vec3 P{rng(-5, 7), rng(-8, 4), rng(-3, 3)};
    Vec3 radial = P - C;
    radial.z = 0;
    if (norm(radial) < 1e-3) continue;
    const Vec3 foot = C + R * (radial / norm(radial));
    const auto r = closest_point(circle, P);
    CHECK(distance(r.point, foot) < 1e-9);
    CHECK(std::abs(r.distance - distance(P, foot)) < 1e-9);
  }
  // centre: every point is a minimiser
  CHECK(std::abs(closest_point(circle, C).distance - R) < 1e-12);
}

TEST_CASE("curve closest point on the domain boundary", "[closest]") {
  const auto line = make_line({0, 0, 0}, {2, 0, 0});
  const auto r = closest_point(line, {5, 1, 0});
  CHECK(r.u == 1.0);
  CHECK(std::abs(r.distance - std::sqrt(10.0)) < 1e-12);
  CHECK(closest_point(line, {-3, 0, 4}).u == 0.0);
}

TEST_CASE("local curve projection", "[closest]") {
  const auto c = wavy_cubic();
  const Vec3 target = c.evaluate(0.62);
  const auto r = closest_point_local(c, target + Vec3{0, 0, 1e-3}, 0.6);
  CHECK(std::abs(r.u - 0.62) < 1e-4);
  CHECK(closest_point_local(c, target, 0.58).distance < 1e-10);
}

// ------------------------------------------------------------------ surfaces

TEST_CASE("closest point on a sphere, including poles and the axis", "[closest][surface]") {
  const double R = 2.0;
  const auto meridian = make_arc({0, 0, 0}, {0, 0, -1}, {1, 0, 0}, R, 0.0, std::numbers::pi);
  const auto sphere = revolve(meridian, {0, 0, 0}, {0, 0, 1}, 2 * std::numbers::pi);
  Rng rng;
  std::vector<Vec3> points{{0, 0, 5}, {0, 0, -0.3}, {0, 0, 1.9}, {1e-9, 0, 3}};
  for (int k = 0; k < 30; ++k) points.push_back({rng(-4, 4), rng(-4, 4), rng(-4, 4)});
  for (const auto& P : points) {
    const double d = norm(P);
    if (d < 1e-3) continue;
    const auto r = closest_point(sphere, P);
    CHECK(std::abs(r.distance - std::abs(d - R)) < 1e-9);
    CHECK(distance(r.point, R * (P / d)) < 1e-6);
  }
}

TEST_CASE("surface point inversion", "[closest][surface]") {
  const auto torus = revolve(make_circle({5, 0, 0}, {0, 1, 0}, 1.5), {0, 0, 0}, {0, 0, 1},
                             2 * std::numbers::pi);
  const auto patch = rational_patch();
  for (const auto* s : {&torus, &patch})
    for (double u = 0.05; u < 1.0; u += 0.13)
      for (double v = 0.03; v < 1.0; v += 0.17) {
        const Vec3 P = s->evaluate(u, v);
        const auto r = closest_point(*s, P);
        CHECK(r.distance < 1e-9);
        CHECK(distance(s->evaluate(r.u, r.v), P) < 1e-9);
      }
}

TEST_CASE("surface closest point is never worse than brute force", "[closest][surface]") {
  const auto patch = rational_patch();
  const auto torus = revolve(make_circle({5, 0, 0}, {0, 1, 0}, 1.5), {0, 0, 0}, {0, 0, 1},
                             3 * std::numbers::pi / 2);
  Rng rng;
  for (int k = 0; k < 12; ++k) {
    const Vec3 P{rng(-1, 5), rng(-1, 4), rng(-3, 3)};
    CHECK(closest_point(patch, P).distance <= brute_force(patch, P) + 1e-9);
    const Vec3 Q{rng(-8, 8), rng(-8, 8), rng(-3, 3)};
    CHECK(closest_point(torus, Q).distance <= brute_force(torus, Q) + 1e-9);
  }
}

TEST_CASE("surface closest point on the boundary and local projection", "[closest][surface]") {
  const auto rect = make_rectangle({0, 0, 0}, {2, 0, 0}, {0, 1, 0});
  const auto r = closest_point(rect, {3, 0.5, 1});
  CHECK(r.u == 1.0);
  CHECK(std::abs(r.v - 0.5) < 1e-12);
  CHECK(std::abs(r.distance - std::sqrt(2.0)) < 1e-12);

  const auto patch = rational_patch();
  const Vec3 P = patch.evaluate(0.45, 0.35);
  const auto l = closest_point_local(patch, P, 0.4, 0.3);
  CHECK(l.distance < 1e-10);
}
