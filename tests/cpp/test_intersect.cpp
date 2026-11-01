#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <numbers>
#include <vector>

#include "nurbforge/analysis/intersect.hpp"
#include "nurbforge/geometry/primitives.hpp"
#include "nurbforge/modeling/construct.hpp"

using namespace nf;
constexpr double pi = std::numbers::pi;

namespace {

NurbsCurve wavy_cubic() {
  return NurbsCurve::from_points(
      3, {{0, 0.5, 0}, {1, 3, 0}, {2, -2, 0}, {3, 3, 0}, {4, -1, 0}, {5, 2, 0}, {6, -0.5, 0}});
}

// Roots of f on [a, b] by sign changes on a fine grid and bisection.
template <class F>
std::vector<double> roots(F f, double a, double b, int n = 20000) {
  std::vector<double> r;
  double x0 = a, f0 = f(a);
  for (int i = 1; i <= n; ++i) {
    const double x1 = a + (b - a) * i / n, f1 = f(x1);
    if (f0 == 0.0) r.push_back(x0);
    if (f0 * f1 < 0.0) {
      double lo = x0, hi = x1;
      for (int k = 0; k < 100; ++k) {
        const double m = 0.5 * (lo + hi);
        (f(m) * f0 <= 0.0 ? hi : lo) = m;
      }
      r.push_back(0.5 * (lo + hi));
    }
    x0 = x1;
    f0 = f1;
  }
  return r;
}

NurbsSurface graph_patch() {
  std::vector<Vec3> P;
  for (int i = 0; i < 6; ++i)
    for (int j = 0; j < 5; ++j) P.push_back({double(i), double(j), std::sin(1.7 * i) * std::cos(1.3 * j)});
  return {3, 3, clamped_uniform_knots(3, 6), clamped_uniform_knots(3, 5), 6, 5, P};
}

}  // namespace

TEST_CASE("line / line intersections", "[intersect]") {
  const auto r = intersect(make_line({0, 0, 0}, {2, 2, 0}), make_line({0, 2, 0}, {2, 0, 0}));
  REQUIRE(r.size() == 1);
  CHECK(distance(r[0].point, {1, 1, 0}) < 1e-14);
  CHECK(std::abs(r[0].s - 0.5) < 1e-14);
  CHECK(intersect(make_line({0, 0, 0}, {1, 0, 0}), make_line({0, 1, 1}, {1, 1, -1})).empty());  // skew
  CHECK(intersect(make_line({0, 0, 0}, {1, 0, 0}), make_line({0, 1, 0}, {1, 1, 0})).empty());   // parallel
  const auto e = intersect(make_line({0, 0, 0}, {1, 0, 0}), make_line({1, 0, 0}, {1, 1, 0}));  // shared end
  REQUIRE(e.size() == 1);
  CHECK(e[0].s == 1.0);
  CHECK(e[0].t == 0.0);
}

TEST_CASE("circle / line and circle / circle", "[intersect]") {
  const auto circle = make_circle({0, 0, 0}, {0, 0, 1}, 2.0);
  const auto r = intersect(circle, make_line({-3, 1, 0}, {3, 1, 0}));
  REQUIRE(r.size() == 2);
  for (const auto& x : r) {
    CHECK(std::abs(x.point.y - 1.0) < 1e-12);
    CHECK(std::abs(std::abs(x.point.x) - std::sqrt(3.0)) < 1e-12);
  }
  const auto t = intersect(circle, make_line({-3, 2, 0}, {3, 2, 0}));  // tangent
  REQUIRE(t.size() == 1);
  CHECK(distance(t[0].point, {0, 2, 0}) < 1e-5);
  CHECK(intersect(circle, make_line({-3, 2.1, 0}, {3, 2.1, 0})).empty());
  const auto cc = intersect(circle, make_circle({2, 0, 0}, {0, 0, 1}, 2.0));
  REQUIRE(cc.size() == 2);
  for (const auto& x : cc) {
    CHECK(std::abs(x.point.x - 1.0) < 1e-12);
    CHECK(std::abs(std::abs(x.point.y) - std::sqrt(3.0)) < 1e-12);
  }
}

TEST_CASE("many intersections of a wavy curve with a line", "[intersect]") {
  const auto c = wavy_cubic();
  for (double y0 : {0.0, 0.7, -0.4}) {
    const auto r = intersect(c, make_line({-1, y0, 0}, {7, y0, 0}));
    const auto ref = roots([&](double t) { return c.evaluate(t).y - y0; }, 0.0, 1.0);
    REQUIRE(r.size() == ref.size());
    for (std::size_t i = 0; i < r.size(); ++i) {
      CHECK(std::abs(r[i].s - ref[i]) < 1e-9);
      CHECK(r[i].distance < 1e-10);
    }
  }
}

TEST_CASE("line / sphere", "[intersect][surface]") {
  const auto meridian = make_arc({0, 0, 0}, {0, 0, -1}, {1, 0, 0}, 2.0, 0.0, pi);
  const auto sphere = revolve(meridian, {0, 0, 0}, {0, 0, 1}, 2 * pi);
  const auto r = intersect(make_line({-5, 0.5, 0.3}, {5, 0.5, 0.3}), sphere);
  REQUIRE(r.size() == 2);
  for (const auto& x : r) CHECK(std::abs(norm(x.point) - 2.0) < 1e-12);
  CHECK(r[0].t < r[1].t);
  // through both poles (on the axis)
  const auto poles = intersect(make_line({0, 0, -3}, {0, 0, 3}), sphere);
  REQUIRE(poles.size() == 2);
  CHECK(distance(poles[0].point, {0, 0, -2}) < 1e-12);
  CHECK(distance(poles[1].point, {0, 0, 2}) < 1e-12);
  const auto tangent = intersect(make_line({-5, 2, 0}, {5, 2, 0}), sphere);
  REQUIRE(tangent.size() == 1);
  CHECK(distance(tangent[0].point, {0, 2, 0}) < 1e-5);
  CHECK(intersect(make_line({-5, 2.2, 0}, {5, 2.2, 0}), sphere).empty());
}

TEST_CASE("vertical lines through a graph surface hit it exactly once", "[intersect][surface]") {
  const auto s = graph_patch();
  for (double x = 0.3; x < 5.0; x += 0.9)
    for (double y = 0.2; y < 4.0; y += 0.7) {
      const auto r = intersect(make_line({x, y, -5}, {x, y, 5}), s);
      REQUIRE(r.size() == 1);
      const Vec3 p = s.evaluate(r[0].u, r[0].v);
      CHECK(std::abs(p.x - x) < 1e-10);
      CHECK(std::abs(p.y - y) < 1e-10);
    }
}

TEST_CASE("curves crossing surfaces several times", "[intersect][surface]") {
  // axis line through a torus lying in the xz plane: four hits
  const auto torus = revolve(make_circle({5, 0, 0}, {0, 1, 0}, 1.5), {0, 0, 0}, {0, 0, 1}, 2 * pi);
  const auto r = intersect(make_line({-8, 0, 0}, {8, 0, 0}), torus);
  REQUIRE(r.size() == 4);
  const double xs[] = {-6.5, -3.5, 3.5, 6.5};
  for (std::size_t i = 0; i < 4; ++i) CHECK(distance(r[i].point, {xs[i], 0, 0}) < 1e-12);

  // tilted circle against a plane: two hits
  const auto circle = make_circle({0, 0, 0}, {1, 0, 1}, 3.0);
  const auto plane = make_rectangle({-5, -5, 0.5}, {10, 0, 0}, {0, 10, 0});
  const auto cp = intersect(circle, plane);
  REQUIRE(cp.size() == 2);
  for (const auto& x : cp) {
    CHECK(std::abs(x.point.z - 0.5) < 1e-12);
    CHECK(std::abs(norm(x.point) - 3.0) < 1e-12);
  }

  // space curve oscillating through the wavy patch: compare with sign changes of z_c - z_s
  const auto s = graph_patch();
  std::vector<Vec3> P;
  for (int k = 0; k <= 12; ++k) P.push_back({0.2 + 0.38 * k, 0.5 + 0.25 * k, (k % 2 == 0 ? 1.5 : -1.5)});
  const auto c = NurbsCurve::from_points(3, P);
  const auto hits = intersect(c, s);
  // reference: the patch is a graph with x(u) and y(v) monotone, so project each curve point
  const auto height_gap = [&](double t) {
    const Vec3 q = c.evaluate(t);
    // find (u, v) with S(u, v).xy = q.xy by bisection-free Newton (graph is nearly linear in u, v)
    double u = q.x / 5.0, v = q.y / 4.0;
    for (int it = 0; it < 50; ++it) {
      const auto d = s.first_derivatives(u, v);
      const double ex = d[0].x - q.x, ey = d[0].y - q.y;
      const double a = d[1].x, b = d[2].x, cc = d[1].y, dd = d[2].y, det = a * dd - b * cc;
      u -= (dd * ex - b * ey) / det;
      v -= (-cc * ex + a * ey) / det;
    }
    return q.z - s.evaluate(u, v).z;
  };
  const auto ref = roots(height_gap, 0.0, 1.0, 4000);
  REQUIRE(hits.size() == ref.size());
  for (std::size_t i = 0; i < hits.size(); ++i) CHECK(std::abs(hits[i].t - ref[i]) < 1e-8);
}
