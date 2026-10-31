#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "nurbforge/core/error.hpp"
#include "nurbforge/nurbs/basis.hpp"

using Catch::Approx;
using namespace nf;

TEST_CASE("find_span (NURBS Book Ex. 2.3)", "[basis]") {
  const KnotVector U{0, 0, 0, 1, 2, 3, 4, 4, 5, 5, 5};
  CHECK(find_span(7, 2, 2.5, U) == 4);
  CHECK(find_span(7, 2, 0.0, U) == 2);
  CHECK(find_span(7, 2, 5.0, U) == 7);
}

TEST_CASE("basis functions form a partition of unity", "[basis]") {
  const KnotVector U{0, 0, 0, 0, 0.2, 0.5, 0.7, 1, 1, 1, 1};
  for (double u = 0.0; u <= 1.0; u += 0.05) {
    const auto s = find_span(6, 3, u, U);
    const auto N = basis_functions(s, u, 3, U);
    double sum = 0;
    for (double n : N) { sum += n; CHECK(n >= -1e-14); }
    CHECK(sum == Approx(1.0));
  }
}

TEST_CASE("basis derivatives sum to zero", "[basis]") {
  const KnotVector U{0, 0, 0, 1, 2, 3, 4, 4, 5, 5, 5};
  const auto s = find_span(7, 2, 2.5, U);
  const auto D = basis_function_derivatives(s, 2.5, 2, 2, U);
  double d1 = 0, d2 = 0;
  for (std::size_t j = 0; j < 3; ++j) { d1 += D[1][j]; d2 += D[2][j]; }
  CHECK(d1 == Approx(0.0).margin(1e-12));
  CHECK(d2 == Approx(0.0).margin(1e-12));
  // values from the NURBS Book, Ex. 2.4: N_{2,2}(5/2)=1/8, N_{3,2}=6/8, N_{4,2}=1/8
  CHECK(D[0][0] == Approx(0.125));
  CHECK(D[0][1] == Approx(0.75));
  CHECK(D[0][2] == Approx(0.125));
}

TEST_CASE("invalid knot vectors are rejected", "[basis]") {
  CHECK_THROWS_AS(validate_knots({0, 0, 1, 0.5, 1}, 1, 3), GeometryError);
  CHECK_THROWS_AS(validate_knots({0, 0, 1, 1}, 1, 3), GeometryError);
  CHECK_NOTHROW(validate_knots(clamped_uniform_knots(3, 7), 3, 7));
}
