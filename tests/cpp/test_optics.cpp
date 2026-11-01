#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <cmath>
#include <cstdint>
#include <numbers>
#include <vector>

#include "nurbforge/core/error.hpp"
#include "nurbforge/optics/surfaces.hpp"

using namespace nf;
constexpr double pi = std::numbers::pi;

namespace {

struct Rng {
  std::uint64_t s = 4242;
  double operator()(double lo, double hi) {
    s = s * 6364136223846793005ULL + 1442695040888963407ULL;
    return lo + (hi - lo) * static_cast<double>(s >> 11) / static_cast<double>(1ULL << 53);
  }
};

// Max |z - sag(r)| over a dense grid of the surface of revolution.
double revolution_error(const NurbsSurface& s, const Asphere& a) {
  double e = 0.0;
  const auto du = s.domain_u(), dv = s.domain_v();
  for (int i = 0; i <= 24; ++i)
    for (int j = 0; j <= 400; ++j) {
      const Vec3 p = s.evaluate(du[0] + (du[1] - du[0]) * i / 24.0, dv[0] + (dv[1] - dv[0]) * j / 400.0);
      e = std::max(e, std::abs(p.z - a.sag(std::hypot(p.x, p.y))));
    }
  return e;
}

}  // namespace

TEST_CASE("asphere sag formula", "[optics]") {
  const Asphere sphere{1.0 / 50.0, 0.0, {}, {}};
  CHECK(std::abs(sphere.sag(20.0) - (50.0 - std::sqrt(2500.0 - 400.0))) < 1e-12);
  const Asphere parabola{1.0 / 40.0, -1.0, {}, {}};
  CHECK(std::abs(parabola.sag(10.0) - 100.0 / 80.0) < 1e-13);
  const Asphere asph{0.0, 0.0, {1e-3, 2e-5}, {0.5, 0.0, 1e-2}};
  CHECK(std::abs(asph.sag(-2.0) - (1e-3 * 16 + 2e-5 * 64 + 0.5 * 2 + 1e-2 * 8)) < 1e-14);
  CHECK(std::abs(sphere.max_radius() - 50.0) < 1e-12);
  CHECK(std::isinf(parabola.max_radius()));
  CHECK_THROWS_AS(sphere.sag(51.0), GeometryError);
}

TEST_CASE("conic surfaces are exact", "[optics]") {
  const double R = 50.0;
  const double k = GENERATE(0.0, -1.0, -2.5, -0.5, 0.6);
  const double r_max = GENERATE(10.0, 30.0);
  for (double c : {1.0 / R, -1.0 / R}) {
    const Asphere a{c, k, {}, {}};
    const auto s = asphere_surface(a, r_max);
    CHECK(asphere_profile(a, r_max).degree() == 2);
    const auto du = s.domain_u(), dv = s.domain_v();
    for (int i = 0; i <= 12; ++i)
      for (int j = 0; j <= 40; ++j) {
        const Vec3 p = s.evaluate(du[0] + (du[1] - du[0]) * i / 12.0, dv[0] + (dv[1] - dv[0]) * j / 40.0);
        const double r2 = p.x * p.x + p.y * p.y;
        // implicit conic: c r^2 - 2 z + (1 + k) c z^2 = 0
        CHECK(std::abs(c * r2 - 2.0 * p.z + (1.0 + k) * c * p.z * p.z) < 1e-13);
      }
    CHECK(std::abs(s.evaluate(0.3, dv[1]).z - a.sag(r_max)) < 1e-12);  // rim
  }
}

TEST_CASE("hemisphere and the conic domain", "[optics]") {
  const Asphere sphere{1.0 / 10.0, 0.0, {}, {}};
  const auto h = asphere_surface(sphere, 10.0);  // vertical tangent at the rim
  const auto dv = h.domain_v();
  for (double v = dv[0]; v <= dv[1]; v += 0.05 * (dv[1] - dv[0])) {
    const Vec3 p = h.evaluate(0.2, v);
    CHECK(std::abs(std::sqrt(p.x * p.x + p.y * p.y + (p.z - 10.0) * (p.z - 10.0)) - 10.0) < 1e-12);
  }
  CHECK_THROWS_AS(asphere_surface(sphere, 10.5), GeometryError);
}

TEST_CASE("even and odd aspheres within tolerance", "[optics]") {
  const Asphere even{1.0 / 40.0, -0.8, {1e-6, -2e-9, 3e-12}, {}};
  for (double tol : {1e-6, 1e-8}) {
    const auto s = asphere_surface(even, 15.0, tol);
    CHECK(revolution_error(s, even) <= 1.5 * tol);
  }
  // the vertex normal is exactly the axis (symmetric fit)
  const auto s = asphere_surface(even, 15.0);
  CHECK(std::abs(std::abs(s.normal(0.3, 0.0).z) - 1.0) < 1e-12);

  const Asphere odd{1.0 / 60.0, 0.0, {}, {0.0, 0.0, 2e-5, -1e-7}};
  CHECK(revolution_error(asphere_surface(odd, 12.0), odd) <= 1.5e-7);
  const Asphere cone_tip{0.0, 0.0, {}, {0.05}};  // linear term: genuine cone at the vertex
  CHECK(revolution_error(asphere_surface(cone_tip, 5.0), cone_tip) <= 1e-12);
}

TEST_CASE("Zernike polynomials: known forms and indexing", "[optics]") {
  const double rho = 0.7, th = 0.4;
  const double r2 = rho * rho;
  // Fringe
  CHECK(std::abs(zernike(1, ZernikeIndexing::Fringe, rho, th) - 1.0) < 1e-15);
  CHECK(std::abs(zernike(2, ZernikeIndexing::Fringe, rho, th) - rho * std::cos(th)) < 1e-15);
  CHECK(std::abs(zernike(4, ZernikeIndexing::Fringe, rho, th) - (2 * r2 - 1)) < 1e-14);
  CHECK(std::abs(zernike(7, ZernikeIndexing::Fringe, rho, th) - (3 * r2 * rho - 2 * rho) * std::cos(th)) < 1e-14);
  CHECK(std::abs(zernike(9, ZernikeIndexing::Fringe, rho, th) - (6 * r2 * r2 - 6 * r2 + 1)) < 1e-14);
  CHECK(std::abs(zernike(6, ZernikeIndexing::Fringe, rho, th) - r2 * std::sin(2 * th)) < 1e-14);
  // Noll (normalised)
  CHECK(std::abs(zernike(2, ZernikeIndexing::Noll, rho, th) - 2 * rho * std::cos(th)) < 1e-14);
  CHECK(std::abs(zernike(3, ZernikeIndexing::Noll, rho, th) - 2 * rho * std::sin(th)) < 1e-14);
  CHECK(std::abs(zernike(4, ZernikeIndexing::Noll, rho, th) - std::sqrt(3.0) * (2 * r2 - 1)) < 1e-14);
  CHECK(std::abs(zernike(11, ZernikeIndexing::Noll, rho, th) - std::sqrt(5.0) * (6 * r2 * r2 - 6 * r2 + 1)) < 1e-13);
  int n = 0, m = 0;
  zernike_nm(9, ZernikeIndexing::Noll, n, m);
  CHECK((n == 3 && m == -3));
  zernike_nm(16, ZernikeIndexing::Fringe, n, m);  // spherical 6th order is Z16 in Fringe
  CHECK((n == 6 && m == 0));
}

TEST_CASE("Noll Zernike polynomials are orthonormal on the unit disk", "[optics]") {
  // polar midpoint rule
  const int nr = 100, nt = 128;
  for (std::size_t i = 1; i <= 15; ++i)
    for (std::size_t j = i; j <= 15; ++j) {
      double s = 0.0;
      for (int a = 0; a < nr; ++a)
        for (int b = 0; b < nt; ++b) {
          const double rho = (a + 0.5) / nr, th = 2 * pi * (b + 0.5) / nt;
          s += zernike(i, ZernikeIndexing::Noll, rho, th) * zernike(j, ZernikeIndexing::Noll, rho, th) * rho;
        }
      s *= (1.0 / nr) * (2 * pi / nt) / pi;
      CHECK(std::abs(s - (i == j ? 1.0 : 0.0)) < 2e-3);
    }
}

TEST_CASE("Zernike surface within tolerance", "[optics]") {
  ZernikeSag z;
  z.base = {1.0 / 80.0, -0.5, {}, {}};
  z.norm_radius = 10.0;
  z.coeffs = {0.0, 0.01, -0.02, 0.05, 0.003, 0.0, 0.004, -0.002, 0.01};
  const double tol = 1e-7;
  const auto s = zernike_surface(z, 0.0, tol);
  Rng rng;
  for (int k = 0; k < 400; ++k) {
    const double x = rng(-10, 10), y = rng(-10, 10);
    const Vec3 p = s.evaluate(x, y);
    CHECK(std::abs(p.x - x) < 1e-12);
    CHECK(std::abs(p.y - y) < 1e-12);
    CHECK(std::abs(p.z - z.sag(x, y)) < 2 * tol);
  }
}

TEST_CASE("XY polynomial surfaces", "[optics]") {
  XYPolynomialSag cubic;
  cubic.norm_radius = 5.0;
  cubic.terms = {{2, 0, 0.1}, {1, 1, -0.05}, {0, 3, 0.02}, {2, 1, 0.03}};  // degree <= 3: exact
  const auto s = xy_polynomial_surface(cubic);
  CHECK(s.num_u() == 5);  // the initial two spans suffice
  for (double x = -5; x <= 5; x += 0.5)
    for (double y = -5; y <= 5; y += 0.5) CHECK(std::abs(s.evaluate(x, y).z - cubic.sag(x, y)) < 1e-13);

  XYPolynomialSag high = cubic;
  high.base = {1.0 / 100.0, 0.0, {}, {}};
  high.terms.push_back({4, 2, 0.01});
  high.terms.push_back({0, 6, -0.004});
  const auto h = xy_polynomial_surface(high, 0.0, 1e-8);
  Rng rng;
  for (int k = 0; k < 300; ++k) {
    const double x = rng(-5, 5), y = rng(-5, 5);
    CHECK(std::abs(h.evaluate(x, y).z - high.sag(x, y)) < 2e-8);
  }
}

TEST_CASE("sag_surface rejects undefined sag", "[optics]") {
  const Asphere sphere{1.0 / 10.0, 0.0, {}, {}};
  // the corners of the square (r = 9 * sqrt 2) lie outside the sphere
  CHECK_THROWS_AS(sag_surface([&](double x, double y) { return sphere.sag(std::hypot(x, y)); }, -9, 9, -9, 9),
                  GeometryError);
}
