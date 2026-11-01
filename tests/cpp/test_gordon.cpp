#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <numbers>
#include <vector>

#include "nurbforge/core/error.hpp"
#include "nurbforge/geometry/primitives.hpp"
#include "nurbforge/geometry/transform.hpp"
#include "nurbforge/modeling/gordon.hpp"

using namespace nf;
constexpr double pi = std::numbers::pi;

namespace {

template <class F>
void grid(int n, F f) {
  for (int i = 0; i <= n; ++i)
    for (int j = 0; j <= n; ++j) f(double(i) / n, double(j) / n);
}

// Non-rational bicubic Bezier patch.
NurbsSurface bicubic() {
  std::vector<Vec3> P;
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j) P.push_back({double(i), double(j), std::sin(1.3 * i + 0.4) * std::cos(0.8 * j) + 0.1 * i * j});
  return {3, 3, {0, 0, 0, 0, 1, 1, 1, 1}, {0, 0, 0, 0, 1, 1, 1, 1}, 4, 4, P};
}

// Non-rational multi-span surface.
NurbsSurface multispan() {
  std::vector<Vec3> P;
  for (int i = 0; i < 7; ++i)
    for (int j = 0; j < 6; ++j) P.push_back({double(i), double(j), std::sin(0.9 * i) * std::cos(1.1 * j)});
  return {3, 2, clamped_uniform_knots(3, 7), clamped_uniform_knots(2, 6), 7, 6, P};
}

}  // namespace

TEST_CASE("Coons patch of four lines is the bilinear patch", "[coons]") {
  const Vec3 p00{0, 0, 0}, p10{2, 0, 1}, p01{0, 3, 0}, p11{2, 3, 4};
  const auto s = coons(make_line(p00, p10), make_line(p01, p11), make_line(p00, p01), make_line(p10, p11));
  const auto b = make_bilinear(p00, p10, p01, p11);
  grid(10, [&](double u, double v) { CHECK(distance(s.evaluate(u, v), b.evaluate(u, v)) < 1e-12); });
}

TEST_CASE("Coons patch interpolates curved planar boundaries", "[coons]") {
  const auto c0 = NurbsCurve::from_points(3, {{0, 0, 0}, {1, -0.5, 0}, {2, 0.4, 0}, {3, 0, 0}});
  const auto c1 = NurbsCurve::from_points(2, {{0, 3, 0}, {1.5, 3.8, 0}, {3, 3, 0}});
  const auto d0 = NurbsCurve::from_points(3, {{0, 0, 0}, {-0.6, 1, 0}, {0.4, 2, 0}, {0, 3, 0}});
  const auto d1 = make_line({3, 0, 0}, {3, 3, 0});
  const auto s = coons(c0, c1, d0, d1);
  grid(20, [&](double u, double v) { CHECK(std::abs(s.evaluate(u, v).z) < 1e-14); });
  for (double t = 0.0; t <= 1.0; t += 0.05) {
    CHECK(distance(s.evaluate(t, 0.0), c0.evaluate(t)) < 1e-12);
    CHECK(distance(s.evaluate(t, 1.0), c1.evaluate(t)) < 1e-12);
    CHECK(distance(s.evaluate(0.0, t), d0.evaluate(t)) < 1e-12);
    CHECK(distance(s.evaluate(1.0, t), d1.evaluate(t)) < 1e-12);
  }
  // scrambled orientations give the same patch
  const auto r = coons(c0, c1.reversed(), d0.reversed(), d1.reversed());
  grid(10, [&](double u, double v) { CHECK(distance(r.evaluate(u, v), s.evaluate(u, v)) < 1e-12); });
}

TEST_CASE("Coons patch with rational arc boundaries", "[coons]") {
  const double R = 2.0;
  const auto c0 = make_arc({0, 0, 0}, {1, 0, 0}, {0, 1, 0}, R, 0.0, pi / 2);
  const auto c1 = c0.transformed(Transform::translation({0, 0, 3}));
  const auto d0 = make_line({R, 0, 0}, {R, 0, 3});
  const auto d1 = make_line({0, R, 0}, {0, R, 3});
  const double tol = 1e-7;
  const auto s = coons(c0, c1, d0, d1, tol);
  grid(20, [&](double u, double v) {
    const Vec3 p = s.evaluate(u, v);
    CHECK(std::abs(std::hypot(p.x, p.y) - R) < 2 * tol);  // a quarter cylinder
  });
}

TEST_CASE("Coons patch rejects an open boundary", "[coons]") {
  CHECK_THROWS_AS(coons(make_line({0, 0, 0}, {1, 0, 0}), make_line({0, 1, 0}, {1, 1, 0}),
                        make_line({0, 0, 0}, {0, 1, 0}), make_line({1.1, 0, 0}, {1, 1, 0})),
                  GeometryError);
}

TEST_CASE("Gordon surface reproduces a bicubic patch from a 4x4 network", "[gordon]") {
  const auto S = bicubic();
  const std::vector<double> us{0, 0.3, 0.65, 1}, vs{0, 0.4, 0.7, 1};
  std::vector<NurbsCurve> uc, vc;
  for (double v : vs) uc.push_back(S.iso_curve_v(v));
  for (double u : us) vc.push_back(S.iso_curve_u(u));
  const auto g = gordon(uc, vc);
  grid(20, [&](double u, double v) { CHECK(distance(g.evaluate(u, v), S.evaluate(u, v)) < 1e-10); });
}

TEST_CASE("Gordon surface interpolates a multi-span network, in any order", "[gordon]") {
  const auto S = multispan();
  const std::vector<double> us{0, 0.22, 0.5, 0.81, 1}, vs{0, 0.35, 1};
  std::vector<NurbsCurve> uc, vc;
  for (double v : vs) uc.push_back(S.iso_curve_v(v));
  for (double u : us) vc.push_back(S.iso_curve_u(u));
  const auto g = gordon(uc, vc);
  for (double t = 0.0; t <= 1.0; t += 0.02) {
    for (double v : vs) CHECK(distance(g.evaluate(t, v), S.evaluate(t, v)) < 1e-10);
    for (double u : us) CHECK(distance(g.evaluate(u, t), S.evaluate(u, t)) < 1e-10);
  }
  std::vector<NurbsCurve> uc2{uc[2], uc[0], uc[1]}, vc2{vc[3], vc[1], vc[4], vc[0], vc[2]};
  const auto h = gordon(uc2, vc2);
  grid(10, [&](double u, double v) { CHECK(distance(h.evaluate(u, v), g.evaluate(u, v)) < 1e-10); });
}

TEST_CASE("Gordon surface rejects inconsistent networks", "[gordon]") {
  const auto S = bicubic();
  std::vector<NurbsCurve> uc{S.iso_curve_v(0.0), S.iso_curve_v(0.5), S.iso_curve_v(1.0)};
  std::vector<NurbsCurve> vc{S.iso_curve_u(0.0), S.iso_curve_u(0.5), S.iso_curve_u(1.0)};
  vc[1] = vc[1].transformed(Transform::translation({0, 0, 0.01}));  // misses the u-curves
  CHECK_THROWS_AS(gordon(uc, vc), GeometryError);
  // reparametrised middle u-curve: same geometry (projective reparametrisation of a Bezier
  // curve, weights w_i * 2^i), so it meets every v-curve – but at other parameters
  const auto c = uc[1];
  std::vector<double> w{1, 2, 4, 8};
  std::vector<NurbsCurve> uc3{uc[0], NurbsCurve(c.degree(), c.knots(), c.control_points(), w), uc[2]};
  std::vector<NurbsCurve> vc3{S.iso_curve_u(0.0), S.iso_curve_u(0.5), S.iso_curve_u(1.0)};
  CHECK_THROWS_AS(gordon(uc3, vc3), GeometryError);
  CHECK_THROWS_AS(gordon({uc[0]}, vc3), GeometryError);
}
