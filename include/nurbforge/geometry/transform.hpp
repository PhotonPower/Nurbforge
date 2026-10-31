#pragma once
/// @file transform.hpp
/// Affine transformations in 3D (translation, rotation, scaling, mirroring).
/// NURBS are invariant under affine maps, so transforming control points is exact.

#include <array>

#include "nurbforge/core/vec.hpp"

namespace nf {

class Transform {
 public:
  /// Identity.
  Transform();

  /// Row-major 3x4 matrix [R | t].
  explicit Transform(const std::array<double, 12>& m) : m_(m) {}

  static Transform translation(const Vec3& d);
  /// Rotation by angle (radians) about the axis through `origin` with direction `axis`.
  static Transform rotation(const Vec3& origin, const Vec3& axis, double angle);
  /// Uniform scaling about `center`.
  static Transform scaling(const Vec3& center, double factor);
  /// Non-uniform scaling about `center`.
  static Transform scaling(const Vec3& center, const Vec3& factors);
  /// Reflection in the plane through `point` with normal `normal`.
  static Transform mirror(const Vec3& point, const Vec3& normal);

  [[nodiscard]] Vec3 apply(const Vec3& p) const;
  [[nodiscard]] Vec3 apply_vector(const Vec3& v) const;  ///< ignores translation

  /// Composition: (a * b).apply(p) == a.apply(b.apply(p)).
  friend Transform operator*(const Transform& a, const Transform& b);

  [[nodiscard]] const std::array<double, 12>& matrix() const noexcept { return m_; }

 private:
  std::array<double, 12> m_;
};

}  // namespace nf
