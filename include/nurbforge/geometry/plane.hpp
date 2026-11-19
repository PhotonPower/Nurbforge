#pragma once
/// @file plane.hpp
/// Infinite plane (sections, splitting, projection, fitting).

#include "nurbforge/core/vec.hpp"

namespace nf {

/// Infinite plane through `point` with normal `normal` (normalised internally).
struct Plane {
  Vec3 point;
  Vec3 normal{0, 0, 1};
  /// Signed distance of p from the plane.
  [[nodiscard]] double distance(const Vec3& p) const;
};

}  // namespace nf
