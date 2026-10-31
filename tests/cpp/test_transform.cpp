#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <numbers>

#include "nurbforge/geometry/primitives.hpp"
#include "nurbforge/geometry/transform.hpp"

using Catch::Approx;
using namespace nf;

TEST_CASE("rotation about z", "[transform]") {
  const auto t = Transform::rotation({0, 0, 0}, {0, 0, 1}, std::numbers::pi / 2);
  const auto p = t.apply({1, 0, 5});
  CHECK(p.x == Approx(0).margin(1e-14)); CHECK(p.y == Approx(1)); CHECK(p.z == Approx(5));
}

TEST_CASE("mirror is an involution", "[transform]") {
  const auto m = Transform::mirror({1, 1, 1}, {1, 2, 3});
  const Vec3 p{3, -2, 7};
  CHECK(distance((m * m).apply(p), p) < 1e-12);
}

TEST_CASE("transformed circle remains exact", "[transform]") {
  const auto c = make_circle({0, 0, 0}, {0, 0, 1}, 1.0)
                     .transformed(Transform::rotation({0, 0, 0}, {1, 0, 0}, 0.7) *
                                  Transform::translation({1, 2, 3}));
  const Vec3 center = Transform::rotation({0, 0, 0}, {1, 0, 0}, 0.7).apply({1, 2, 3});
  for (double u = 0; u <= 1; u += 0.1) CHECK(distance(c.evaluate(u), center) == Approx(1.0));
}
