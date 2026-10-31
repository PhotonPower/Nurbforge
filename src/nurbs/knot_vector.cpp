#include "nurbforge/nurbs/knot_vector.hpp"

#include <cmath>
#include <string>

#include "nurbforge/core/error.hpp"
#include "nurbforge/core/vec.hpp"

namespace nf {

void validate_knots(const KnotVector& knots, std::size_t degree, std::size_t n_ctrl) {
  if (degree < 1) throw GeometryError("degree must be >= 1");
  if (n_ctrl < degree + 1)
    throw GeometryError("need at least degree+1 control points (got " + std::to_string(n_ctrl) + ")");
  if (knots.size() != n_ctrl + degree + 1)
    throw GeometryError("knot vector size must be n_ctrl + degree + 1 = " +
                        std::to_string(n_ctrl + degree + 1) + " (got " +
                        std::to_string(knots.size()) + ")");
  for (std::size_t i = 1; i < knots.size(); ++i)
    if (knots[i] < knots[i - 1]) throw GeometryError("knot vector must be non-decreasing");
  if (knots[degree] >= knots[n_ctrl])
    throw GeometryError("knot vector has an empty parameter domain");
  // interior multiplicity must not exceed degree (curve would be disconnected)
  std::size_t mult = 1;
  for (std::size_t i = degree + 1; i < n_ctrl; ++i) {
    mult = (std::abs(knots[i] - knots[i - 1]) < tol::kParametric) ? mult + 1 : 1;
    if (mult > degree) throw GeometryError("interior knot multiplicity exceeds degree");
  }
}

KnotVector clamped_uniform_knots(std::size_t degree, std::size_t n_ctrl) {
  if (n_ctrl < degree + 1) throw GeometryError("need at least degree+1 control points");
  KnotVector k(n_ctrl + degree + 1, 0.0);
  const std::size_t n_inner = n_ctrl - degree - 1;
  for (std::size_t i = 0; i < n_inner; ++i)
    k[degree + 1 + i] = static_cast<double>(i + 1) / static_cast<double>(n_inner + 1);
  for (std::size_t i = n_ctrl; i < k.size(); ++i) k[i] = 1.0;
  return k;
}

std::size_t find_span(std::size_t n, std::size_t degree, double u, const KnotVector& knots) {
  if (u >= knots[n + 1]) return n;
  if (u <= knots[degree]) return degree;
  std::size_t low = degree, high = n + 1;
  std::size_t mid = (low + high) / 2;
  while (u < knots[mid] || u >= knots[mid + 1]) {
    if (u < knots[mid]) high = mid; else low = mid;
    mid = (low + high) / 2;
  }
  return mid;
}

std::size_t knot_multiplicity(const KnotVector& knots, double u) {
  std::size_t m = 0;
  for (double k : knots) if (std::abs(k - u) < tol::kParametric) ++m;
  return m;
}

std::pair<std::vector<double>, std::vector<int>> unique_knots(const KnotVector& knots) {
  std::vector<double> values;
  std::vector<int> mults;
  for (double k : knots) {
    if (!values.empty() && std::abs(k - values.back()) < tol::kParametric) ++mults.back();
    else { values.push_back(k); mults.push_back(1); }
  }
  return {values, mults};
}

}  // namespace nf
