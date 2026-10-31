#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <numbers>
#include <vector>

#include "nurbforge/core/error.hpp"
#include "nurbforge/geometry/primitives.hpp"
#include "nurbforge/modeling/construct.hpp"
#include "nurbforge/nurbs/surface.hpp"

using namespace nf;

namespace {

// Rational bicubic-by-quadratic patch with non-uniform knots and varying weights.
NurbsSurface rational_patch() {
  const std::size_t nu = 5, nv = 4;
  std::vector<Vec3> P;
  std::vector<double> W;
  for (std::size_t i = 0; i < nu; ++i)
    for (std::size_t j = 0; j < nv; ++j) {
      const double x = static_cast<double>(i), y = static_cast<double>(j);
      P.push_back({x, y, std::sin(x) * std::cos(1.3 * y) + 0.1 * x * y});
      W.push_back(1.0 + 0.4 * std::sin(1.7 * x + 0.9 * y) * std::sin(1.7 * x + 0.9 * y));
    }
  return {3, 2, {0, 0, 0, 0, 0.4, 1, 1, 1, 1}, {0, 0, 0, 0.6, 1, 1, 1}, nu, nv, P, W};
}

NurbsSurface sphere(double R) {
  const auto meridian = make_arc({0, 0, 0}, {0, 0, -1}, {1, 0, 0}, R, 0.0, std::numbers::pi);
  return revolve(meridian, {0, 0, 0}, {0, 0, 1}, 2 * std::numbers::pi);
}

// Checks D[k][l] against central differences of the analytic D[k-1][l] (in u) and D[k][l-1]
// (in v) for all 1 <= k + l <= order.
void check_against_finite_differences(const NurbsSurface& s, double u, double v, std::size_t order) {
  const double h = 1e-5;
  const auto D = s.derivatives(u, v, order);
  const auto Dup = s.derivatives(u + h, v, order), Dum = s.derivatives(u - h, v, order);
  const auto Dvp = s.derivatives(u, v + h, order), Dvm = s.derivatives(u, v - h, order);
  for (std::size_t k = 0; k <= order; ++k)
    for (std::size_t l = 0; k + l <= order; ++l) {
      if (k + l == 0) continue;
      const double scale = 1.0 + norm(D[k][l]);
      if (k > 0) {
        const Vec3 fd = (Dup[k - 1][l] - Dum[k - 1][l]) / (2 * h);
        CHECK(distance(fd, D[k][l]) < 1e-5 * scale);
      }
      if (l > 0) {
        const Vec3 fd = (Dvp[k][l - 1] - Dvm[k][l - 1]) / (2 * h);
        CHECK(distance(fd, D[k][l]) < 1e-5 * scale);
      }
    }
}

}  // namespace

TEST_CASE("surface derivatives agree with finite differences", "[derivatives]") {
  const auto s = rational_patch();
  for (double u : {0.13, 0.55, 0.87})  // spans on both sides of the interior knots
    for (double v : {0.21, 0.74}) check_against_finite_differences(s, u, v, 4);
  const auto torus = revolve(make_circle({5, 0, 0}, {0, 1, 0}, 1.5), {0, 0, 0}, {0, 0, 1},
                             2 * std::numbers::pi);
  check_against_finite_differences(torus, 0.33, 0.61, 3);
}

TEST_CASE("derivatives: shape of the result and consistency", "[derivatives]") {
  const auto s = rational_patch();
  const auto D = s.derivatives(0.3, 0.7, 2);
  REQUIRE(D.size() == 3);
  REQUIRE(D[0].size() == 3);
  CHECK(distance(D[0][0], s.evaluate(0.3, 0.7)) < 1e-14);
  const auto f = s.first_derivatives(0.3, 0.7);
  CHECK(distance(f[1], D[1][0]) < 1e-14);
  CHECK(distance(f[2], D[0][1]) < 1e-14);
  CHECK(norm(D[2][1]) == 0.0);  // k + l > order stays zero
  CHECK(norm(D[2][2]) == 0.0);
}

TEST_CASE("bilinear patch: exact second derivatives, higher orders vanish", "[derivatives]") {
  const Vec3 p00{0, 0, 0}, p10{2, 0, 1}, p01{0, 3, 0}, p11{2, 3, 4};
  const auto s = make_bilinear(p00, p10, p01, p11);
  const auto D = s.derivatives(0.4, 0.8, 3);
  CHECK(norm(D[2][0]) < 1e-14);
  CHECK(norm(D[0][2]) < 1e-14);
  CHECK(distance(D[1][1], p00 - p10 - p01 + p11) < 1e-14);
  CHECK(norm(D[2][1]) < 1e-14);
  CHECK(norm(D[3][0]) < 1e-14);
}

TEST_CASE("sphere: normals are radial, including at the poles", "[derivatives][normal]") {
  const double R = 2.5;
  const auto s = sphere(R);
  // orientation of the parametrisation, taken at a regular point
  const double orient = dot(s.normal(0.3, 0.5), s.evaluate(0.3, 0.5)) > 0 ? 1.0 : -1.0;
  for (double u : {0.0, 0.17, 0.5, 0.83}) {
    for (double v : {0.0, 1e-12, 0.3, 0.7, 1.0 - 1e-12, 1.0}) {
      const Vec3 n = s.normal(u, v);
      const Vec3 radial = s.evaluate(u, v) / R;
      CHECK(distance(n, orient * radial) < 1e-6);
    }
  }
}

TEST_CASE("pole normal of a paraboloid-like cap", "[normal]") {
  // quarter circle from the axis outwards, revolved: smooth cap with a pole at v = 0
  const auto profile = make_arc({0, 0, 1}, {0, 0, 1}, {1, 0, 0}, 1.0, 0.0, std::numbers::pi / 2);
  const auto cap = revolve(profile, {0, 0, 0}, {0, 0, 1}, 2 * std::numbers::pi);
  const Vec3 near = cap.normal(0.2, 1e-4);
  for (double u : {0.0, 0.2, 0.6}) {
    const Vec3 n = cap.normal(u, 0.0);
    CHECK(std::abs(std::abs(n.z) - 1.0) < 1e-9);
    CHECK(dot(n, near) > 0.99);  // same orientation as the neighbourhood
  }
}

TEST_CASE("normal is undefined where both partials vanish", "[normal]") {
  const auto s = make_bilinear({0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {1, 1, 1});
  CHECK_THROWS_AS(s.normal(0.0, 0.0), GeometryError);
}
