#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstdint>
#include <numbers>
#include <vector>

#include "nurbforge/core/error.hpp"
#include "nurbforge/fitting/fit.hpp"
#include "nurbforge/geometry/primitives.hpp"
#include "nurbforge/modeling/construct.hpp"

using namespace nf;
constexpr double pi = std::numbers::pi;

namespace {

struct Rng {
  std::uint64_t s{12345};
  double operator()() {
    s = s * 6364136223846793005ULL + 1442695040888963407ULL;
    return static_cast<double>(s >> 11) / 9007199254740992.0;
  }
};

std::vector<Vec3> height_field(std::size_t n, double noise, Rng& rng) {
  std::vector<Vec3> p;
  for (std::size_t k = 0; k < n; ++k) {
    const double x = 4.0 * rng(), y = 3.0 * rng();
    p.push_back({x, y, 0.2 * std::sin(x) * std::cos(y) + noise * (rng() - 0.5)});
  }
  // corners so that the bounding rectangle is the full domain
  for (const auto& [x, y] : {std::pair{0.0, 0.0}, {4.0, 0.0}, {0.0, 3.0}, {4.0, 3.0}})
    p.push_back({x, y, 0.2 * std::sin(x) * std::cos(y)});
  return p;
}

}  // namespace

TEST_CASE("fit surface to scattered points", "[fitting][surface]") {
  Rng rng;
  SECTION("height field: smooth data are reproduced closely") {
    const auto pts = height_field(2000, 0.0, rng);
    const auto f = fit_surface(pts);
    CHECK(f.surface.num_u() == 22);
    CHECK(f.max_error <= 1e-4);
    CHECK(f.parameters.size() == pts.size());
    // the surface follows the function between the data points too
    double dev = 0.0;
    for (int i = 2; i <= 18; ++i)
      for (int j = 2; j <= 18; ++j) {
        const Vec3 q = f.surface.evaluate(i / 20.0, j / 20.0);
        dev = std::max(dev, std::abs(q.z - 0.2 * std::sin(q.x) * std::cos(q.y)));
      }
    CHECK(dev <= 1e-4);
  }
  SECTION("points on a plane: exact up to rounding, any smoothing") {
    std::vector<Vec3> pts;
    for (int k = 0; k < 300; ++k) {
      const double a = 2 * rng() - 1, b = 2 * rng() - 1;
      pts.push_back(Vec3{1, 2, 3} + a * Vec3{1, 1, 0} + b * Vec3{0, 0.5, 2});
    }
    for (const double lambda : {0.0, 1e-4, 1.0}) {
      SurfaceFitOptions o;
      o.smoothing = lambda;
      o.num_u = o.num_v = 6;
      CHECK(fit_surface(pts, o).max_error <= 1e-10);
    }
  }
  SECTION("noise is smoothed out") {
    const auto noisy = height_field(3000, 0.02, rng);
    SurfaceFitOptions loose, tight;
    loose.smoothing = 1.0;
    tight.smoothing = 0.0;
    const auto a = fit_surface(noisy, loose), b = fit_surface(noisy, tight);
    // distance from the noise-free function at the data parameters' xy positions
    const auto truth = [](const Vec3& p) { return 0.2 * std::sin(p.x) * std::cos(p.y); };
    double ea = 0.0, eb = 0.0;
    for (int i = 2; i <= 18; ++i)
      for (int j = 2; j <= 18; ++j) {
        const Vec3 pa = a.surface.evaluate(i / 20.0, j / 20.0);
        const Vec3 pb = b.surface.evaluate(i / 20.0, j / 20.0);
        ea = std::max(ea, std::abs(pa.z - truth(pa)));
        eb = std::max(eb, std::abs(pb.z - truth(pb)));
      }
    CHECK(ea < eb);
    CHECK(ea <= 0.01);
  }
  SECTION("boundary curves are interpolated exactly") {
    const auto f0 = [](double x, double y) { return 0.2 * std::sin(pi * x) * std::sin(pi * y); };
    std::vector<Vec3> pts;
    for (int k = 0; k < 1500; ++k) {
      const double x = rng(), y = rng();
      pts.push_back({x, y, f0(x, y)});
    }
    SurfaceFitOptions o;
    o.boundary = {make_line({0, 0, 0}, {1, 0, 0}), make_line({0, 1, 0}, {1, 1, 0}),
                  make_line({0, 0, 0}, {0, 1, 0}), make_line({1, 0, 0}, {1, 1, 0})};
    const auto f = fit_surface(pts, o);
    CHECK(f.max_error <= 1e-4);
    for (int k = 0; k <= 10; ++k) {
      const double t = k / 10.0;
      CHECK(std::abs(f.surface.evaluate(t, 0.0).z) <= 1e-15);
      CHECK(std::abs(f.surface.evaluate(1.0, t).z) <= 1e-15);
    }
  }
  SECTION("base surface for data that are not a height field") {
    // half cylinder r = 1 around the x axis, points with z from -1 to 1 (folds over in z)
    const auto base = revolve(make_line({0, 0, -1}, {3, 0, -1}), {0, 0, 0}, {1, 0, 0}, pi);
    std::vector<Vec3> pts;
    for (int k = 0; k < 1500; ++k) {
      const double x = 3 * rng(), a = pi * rng();
      const double r = 1.0 + 0.05 * std::sin(2 * x) * std::sin(2 * a);
      pts.push_back({x, r * std::sin(a), -r * std::cos(a)});
    }
    SurfaceFitOptions o;
    o.base = base;
    const auto f = fit_surface(pts, o);
    CHECK(f.max_error <= 1e-3);
    CHECK(f.rms_error <= f.max_error);
  }
  SECTION("errors") {
    CHECK_THROWS_AS(fit_surface({{0, 0, 0}, {1, 0, 0}}), GeometryError);
    CHECK_THROWS_AS(fit_surface({{0, 0, 0}, {1, 0, 0}, {2, 0, 0}, {3, 0, 0}}), GeometryError);
    SurfaceFitOptions o;
    o.boundary = {make_line({0, 0, 0}, {1, 0, 0})};
    CHECK_THROWS_AS(fit_surface(height_field(10, 0, rng), o), GeometryError);
  }
}
