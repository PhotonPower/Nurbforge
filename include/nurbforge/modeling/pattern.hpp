#pragma once
/// @file pattern.hpp
/// Array (pattern) functions: linear, rectangular grid and polar copies.

#include <vector>

#include "nurbforge/geometry/transform.hpp"
#include "nurbforge/nurbs/curve.hpp"
#include "nurbforge/nurbs/surface.hpp"

namespace nf {

/// Transforms for `count` copies along `step` (first copy = identity).
[[nodiscard]] std::vector<Transform> linear_pattern(const Vec3& step, std::size_t count);

/// Rectangular grid: count_a x count_b copies (row-major).
[[nodiscard]] std::vector<Transform> grid_pattern(const Vec3& step_a, std::size_t count_a,
                                                  const Vec3& step_b, std::size_t count_b);

/// `count` copies rotated about an axis, evenly spread over `total_angle` (radians).
/// For a full circle (2*pi) the last copy does not coincide with the first.
[[nodiscard]] std::vector<Transform> polar_pattern(const Vec3& axis_point, const Vec3& axis_dir,
                                                   std::size_t count, double total_angle);

/// Applies a set of transforms to a geometry object.
template <class Geometry>
[[nodiscard]] std::vector<Geometry> apply_pattern(const Geometry& g, const std::vector<Transform>& ts) {
  std::vector<Geometry> out;
  out.reserve(ts.size());
  for (const auto& t : ts) out.push_back(g.transformed(t));
  return out;
}

}  // namespace nf
