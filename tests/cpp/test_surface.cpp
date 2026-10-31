#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "nurbforge/geometry/primitives.hpp"

using Catch::Approx;
using namespace nf;

TEST_CASE("rectangle evaluation and normal", "[surface]") {
  const auto s = make_rectangle({0, 0, 0}, {2, 0, 0}, {0, 3, 0});
  const auto p = s.evaluate(0.5, 0.5);
  CHECK(p.x == Approx(1.0)); CHECK(p.y == Approx(1.5));
  const auto n = s.normal(0.3, 0.7);
  CHECK(n.z == Approx(1.0));
}

TEST_CASE("iso curves and transpose", "[surface]") {
  const auto s = make_bilinear({0, 0, 0}, {1, 0, 1}, {0, 1, 0}, {1, 1, 2});
  const auto cu = s.iso_curve_u(0.25);
  const auto cv = s.iso_curve_v(0.75);
  CHECK(distance(cu.evaluate(0.75), s.evaluate(0.25, 0.75)) < 1e-12);
  CHECK(distance(cv.evaluate(0.25), s.evaluate(0.25, 0.75)) < 1e-12);
  CHECK(distance(s.transposed().evaluate(0.75, 0.25), s.evaluate(0.25, 0.75)) < 1e-12);
}
