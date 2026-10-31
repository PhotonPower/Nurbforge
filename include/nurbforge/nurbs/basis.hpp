#pragma once
/// @file basis.hpp
/// B-spline basis functions (Piegl & Tiller, "The NURBS Book", 2nd ed.).

#include <cstddef>
#include <vector>

#include "nurbforge/nurbs/knot_vector.hpp"

namespace nf {

/// Non-vanishing basis functions N_{span-p..span, p}(u)  (A2.2). Size p+1.
[[nodiscard]] std::vector<double> basis_functions(std::size_t span, double u, std::size_t degree,
                                                  const KnotVector& knots);

/// Basis functions and their derivatives up to order n  (A2.3).
/// Result[k][j] = k-th derivative of N_{span-p+j, p}(u). Size (n+1) x (p+1).
[[nodiscard]] std::vector<std::vector<double>> basis_function_derivatives(
    std::size_t span, double u, std::size_t degree, std::size_t n, const KnotVector& knots);

}  // namespace nf
