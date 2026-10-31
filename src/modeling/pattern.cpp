#include "nurbforge/modeling/pattern.hpp"

#include <cmath>
#include <numbers>

namespace nf {

std::vector<Transform> linear_pattern(const Vec3& step, std::size_t count) {
  std::vector<Transform> out;
  out.reserve(count);
  for (std::size_t i = 0; i < count; ++i)
    out.push_back(Transform::translation(static_cast<double>(i) * step));
  return out;
}

std::vector<Transform> grid_pattern(const Vec3& step_a, std::size_t count_a, const Vec3& step_b,
                                    std::size_t count_b) {
  std::vector<Transform> out;
  out.reserve(count_a * count_b);
  for (std::size_t i = 0; i < count_a; ++i)
    for (std::size_t j = 0; j < count_b; ++j)
      out.push_back(Transform::translation(static_cast<double>(i) * step_a +
                                           static_cast<double>(j) * step_b));
  return out;
}

std::vector<Transform> polar_pattern(const Vec3& axis_point, const Vec3& axis_dir,
                                     std::size_t count, double total_angle) {
  std::vector<Transform> out;
  if (count == 0) return out;
  const bool full = std::abs(total_angle - 2.0 * std::numbers::pi) < tol::kAngular;
  const double denom = static_cast<double>(full ? count : (count > 1 ? count - 1 : 1));
  const double step = total_angle / denom;
  for (std::size_t i = 0; i < count; ++i)
    out.push_back(Transform::rotation(axis_point, axis_dir, step * static_cast<double>(i)));
  return out;
}

}  // namespace nf
