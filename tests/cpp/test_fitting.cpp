#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <algorithm>
#include <cmath>
#include <numbers>
#include <vector>

#include "nurbforge/core/error.hpp"
#include "nurbforge/fitting/fit.hpp"
#include "nurbforge/nurbs/basis.hpp"

using namespace nf;
constexpr double pi = std::numbers::pi;

namespace {

std::vector<Vec3> wiggly_points(std::size_t n) {
  std::vector<Vec3> Q;
  for (std::size_t k = 0; k < n; ++k) {
    const double t = static_cast<double>(k);
    Q.push_back({t, std::sin(1.3 * t) + 0.2 * t, std::cos(0.7 * t)});
  }
  return Q;
}

double max_residual(const NurbsCurve& c, const std::vector<Vec3>& Q, const std::vector<double>& ub) {
  double r = 0.0;
  for (std::size_t k = 0; k < Q.size(); ++k) r = std::max(r, distance(c.evaluate(ub[k]), Q[k]));
  return r;
}

}  // namespace

TEST_CASE("data parameters", "[fitting]") {
  const std::vector<Vec3> Q{{0, 0, 0}, {1, 0, 0}, {3, 0, 0}, {6, 0, 0}};
  const auto u = data_parameters(Q, Parametrization::Uniform);
  const auto c = data_parameters(Q, Parametrization::ChordLength);
  const auto s = data_parameters(Q, Parametrization::Centripetal);
  CHECK(u == std::vector<double>{0.0, 1.0 / 3, 2.0 / 3, 1.0});
  CHECK(std::abs(c[1] - 1.0 / 6) < 1e-15);
  CHECK(std::abs(c[2] - 0.5) < 1e-15);
  const double tot = 1 + std::sqrt(2.0) + std::sqrt(3.0);
  CHECK(std::abs(s[2] - (1 + std::sqrt(2.0)) / tot) < 1e-15);
  CHECK_THROWS_AS(data_parameters({{0, 0, 0}, {0, 0, 0}, {1, 0, 0}}, Parametrization::ChordLength), GeometryError);
  CHECK_THROWS_AS(data_parameters({{1, 1, 1}}, Parametrization::Uniform), GeometryError);
}

TEST_CASE("curve interpolation passes through the points", "[fitting]") {
  const std::size_t p = GENERATE(1, 2, 3, 5);
  const auto method = GENERATE(Parametrization::Uniform, Parametrization::ChordLength, Parametrization::Centripetal);
  const auto Q = wiggly_points(12);
  const auto c = interpolate_curve(Q, p, method);
  CHECK(c.degree() == p);
  CHECK(c.num_control_points() == Q.size());
  CHECK(max_residual(c, Q, data_parameters(Q, method)) < 1e-12);
}

TEST_CASE("curve interpolation reproduces polynomials in the spline space", "[fitting]") {
  std::vector<Vec3> Q;
  for (int k = 0; k <= 10; ++k) {
    const double t = k / 10.0;
    Q.push_back({t, t * t * t, 1 - t * t});
  }
  const auto c = interpolate_curve(Q, 3, Parametrization::Uniform);
  for (double t = 0.0; t <= 1.0; t += 0.01)
    CHECK(distance(c.evaluate(t), Vec3{t, t * t * t, 1 - t * t}) < 1e-12);
}

TEST_CASE("curve interpolation with end derivatives", "[fitting]") {
  const auto Q = wiggly_points(8);
  const Vec3 d0{3, 0, 1}, d1{-2, 4, 0};
  const auto c = interpolate_curve(Q, 3, Parametrization::Centripetal, std::pair{d0, d1});
  CHECK(c.num_control_points() == Q.size() + 2);
  CHECK(max_residual(c, Q, data_parameters(Q, Parametrization::Centripetal)) < 1e-12);
  CHECK(distance(c.derivatives(0.0, 1)[1], d0) < 1e-11);
  CHECK(distance(c.derivatives(1.0, 1)[1], d1) < 1e-11);

  // two points + derivatives, cubic: the Hermite curve
  const Vec3 a{0, 0, 0}, b{1, 1, 0};
  const auto h = interpolate_curve({a, b}, 3, Parametrization::Uniform, std::pair{d0, d1});
  for (double t = 0.0; t <= 1.0; t += 0.05) {
    const double h00 = 2 * t * t * t - 3 * t * t + 1, h10 = t * t * t - 2 * t * t + t;
    const double h01 = -2 * t * t * t + 3 * t * t, h11 = t * t * t - t * t;
    const Vec3 ref = h00 * a + h10 * d0 + h01 * b + h11 * d1;
    CHECK(distance(h.evaluate(t), ref) < 1e-12);
  }
  CHECK_THROWS_AS(interpolate_curve(Q, 1, Parametrization::Uniform, std::pair{d0, d1}), GeometryError);
}

TEST_CASE("curve interpolation input checks", "[fitting]") {
  CHECK_THROWS_AS(interpolate_curve({{0, 0, 0}, {1, 0, 0}, {2, 1, 0}}, 3), GeometryError);
  CHECK_THROWS_AS(interpolate_curve({{0, 0, 0}, {1, 0, 0}, {1, 0, 0}, {2, 1, 0}}, 2), GeometryError);
}

TEST_CASE("least-squares approximation is optimal and keeps the end points", "[fitting]") {
  const auto Q = wiggly_points(40);
  const auto method = Parametrization::ChordLength;
  const auto ub = data_parameters(Q, method);
  const std::size_t p = 3, n_ctrl = 9;
  const auto c = approximate_curve(Q, p, n_ctrl, method);
  CHECK(c.num_control_points() == n_ctrl);
  CHECK(distance(c.control_points().front(), Q.front()) == 0.0);
  CHECK(distance(c.control_points().back(), Q.back()) == 0.0);
  // normal equations: sum_k N_i(u_k) (C(u_k) - Q_k) = 0 for every interior control point
  for (std::size_t i = 1; i + 1 < n_ctrl; ++i) {
    Vec3 g;
    double scale = 0.0;
    for (std::size_t k = 0; k < Q.size(); ++k) {
      const std::size_t span = find_span(n_ctrl - 1, p, ub[k], c.knots());
      if (i + p < span || i > span) continue;
      const double Ni = basis_functions(span, ub[k], p, c.knots())[i + p - span];
      g += Ni * (c.evaluate(ub[k]) - Q[k]);
      scale += Ni * norm(Q[k]);
    }
    CHECK(norm(g) < 1e-12 * (1.0 + scale));
  }
  CHECK_THROWS_AS(approximate_curve(Q, 3, 3), GeometryError);
  CHECK_THROWS_AS(approximate_curve(Q, 3, 41), GeometryError);
}

TEST_CASE("approximation within a tolerance uses few control points", "[fitting]") {
  std::vector<Vec3> Q;
  for (int k = 0; k < 200; ++k) {
    const double a = 1.5 * pi * k / 199.0;
    Q.push_back({std::cos(a), std::sin(a), 0.1 * a});  // helix piece
  }
  for (double tol : {1e-2, 1e-4, 1e-6}) {
    const auto c = approximate_curve_tol(Q, 3, tol);
    const auto ub = data_parameters(Q, Parametrization::Centripetal);
    CHECK(max_residual(c, Q, ub) <= tol);
    CHECK(c.num_control_points() < Q.size());
    const std::size_t n = c.num_control_points();
    if (n > 4) CHECK(max_residual(approximate_curve(Q, 3, n - 1), Q, ub) > tol);
  }
}

TEST_CASE("surface interpolation", "[fitting][surface]") {
  const std::size_t nu = 9, nv = 7;
  std::vector<Vec3> Q;
  for (std::size_t i = 0; i < nu; ++i)
    for (std::size_t j = 0; j < nv; ++j) {
      const double x = static_cast<double>(i) / (nu - 1), y = static_cast<double>(j) / (nv - 1);
      Q.push_back({x, y, x * x * x * y * y});
    }
  // uniform parameters equal (x, y): bicubic x biquadratic data is reproduced exactly
  const auto s = interpolate_surface(Q, nu, nv, 3, 2, Parametrization::Uniform);
  for (double x = 0.0; x <= 1.0; x += 0.1)
    for (double y = 0.0; y <= 1.0; y += 0.1)
      CHECK(distance(s.evaluate(x, y), Vec3{x, y, x * x * x * y * y}) < 1e-12);

  const auto c = interpolate_surface(Q, nu, nv, 3, 3, Parametrization::Centripetal);
  CHECK(c.num_u() == nu);
  CHECK(c.num_v() == nv);
  CHECK_THROWS_AS(interpolate_surface(Q, nu, nv + 1, 3, 3), GeometryError);
}

TEST_CASE("surface interpolation of a sphere grid with pole rows", "[fitting][surface]") {
  const std::size_t nu = 11, nv = 13;  // u: latitude incl. both poles, v: longitude
  std::vector<Vec3> Q;
  for (std::size_t i = 0; i < nu; ++i)
    for (std::size_t j = 0; j < nv; ++j) {
      const double th = pi * static_cast<double>(i) / (nu - 1);
      const double ph = 2 * pi * static_cast<double>(j) / (nv - 1);
      Q.push_back({std::sin(th) * std::cos(ph), std::sin(th) * std::sin(ph), -std::cos(th)});
    }
  const auto s = interpolate_surface(Q, nu, nv, 3, 3, Parametrization::ChordLength);
  // passes through every grid point at the averaged parameters; check a few via the corners/poles
  CHECK(distance(s.evaluate(0.0, 0.3), {0, 0, -1}) < 1e-12);
  CHECK(distance(s.evaluate(1.0, 0.7), {0, 0, 1}) < 1e-12);
  // and stays close to the unit sphere in between
  for (double u = 0.05; u < 1.0; u += 0.1)
    for (double v = 0.0; v <= 1.0; v += 0.1) CHECK(std::abs(norm(s.evaluate(u, v)) - 1.0) < 2e-3);
}

TEST_CASE("surface approximation", "[fitting][surface]") {
  const std::size_t nu = 25, nv = 20;
  std::vector<Vec3> Q;
  for (std::size_t i = 0; i < nu; ++i)
    for (std::size_t j = 0; j < nv; ++j) {
      const double x = static_cast<double>(i) / (nu - 1), y = static_cast<double>(j) / (nv - 1);
      Q.push_back({x, y, 1 + 2 * x - x * x + 3 * x * y * y});
    }
  // data inside the spline space (degree 2 x 2, Bezier) is reproduced exactly
  const auto s = approximate_surface(Q, nu, nv, 2, 2, 3, 3, Parametrization::Uniform);
  for (double x = 0.0; x <= 1.0; x += 0.125)
    for (double y = 0.0; y <= 1.0; y += 0.125)
      CHECK(std::abs(s.evaluate(x, y).z - (1 + 2 * x - x * x + 3 * x * y * y)) < 1e-12);

  const auto w = approximate_surface(Q, nu, nv, 3, 3, 8, 6);
  CHECK(w.num_u() == 8);
  CHECK(w.num_v() == 6);
  CHECK(distance(w.control_point(0, 0), Q.front()) < 1e-15);
  CHECK(distance(w.control_point(7, 5), Q.back()) < 1e-15);
  CHECK_THROWS_AS(approximate_surface(Q, nu, nv, 3, 3, 3, 6), GeometryError);
}
