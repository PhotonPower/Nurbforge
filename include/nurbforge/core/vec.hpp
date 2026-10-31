#pragma once
/// @file vec.hpp
/// Minimal, dependency-free 3D vector type used throughout the kernel.

#include <array>
#include <cmath>
#include <ostream>

#include "nurbforge/core/tolerance.hpp"

namespace nf {

/// Legacy alias for the model-space length tolerance; prefer `nf::tol::kLength`.
inline constexpr double kTolerance = tol::kLength;

struct Vec3 {
  double x{0.0}, y{0.0}, z{0.0};

  constexpr Vec3() = default;
  constexpr Vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}
  explicit constexpr Vec3(const std::array<double, 3>& a) : x(a[0]), y(a[1]), z(a[2]) {}

  [[nodiscard]] constexpr std::array<double, 3> to_array() const { return {x, y, z}; }

  constexpr Vec3& operator+=(const Vec3& o) { x += o.x; y += o.y; z += o.z; return *this; }
  constexpr Vec3& operator-=(const Vec3& o) { x -= o.x; y -= o.y; z -= o.z; return *this; }
  constexpr Vec3& operator*=(double s) { x *= s; y *= s; z *= s; return *this; }
  constexpr Vec3& operator/=(double s) { x /= s; y /= s; z /= s; return *this; }

  friend constexpr Vec3 operator+(Vec3 a, const Vec3& b) { return a += b; }
  friend constexpr Vec3 operator-(Vec3 a, const Vec3& b) { return a -= b; }
  friend constexpr Vec3 operator-(const Vec3& a) { return {-a.x, -a.y, -a.z}; }
  friend constexpr Vec3 operator*(Vec3 a, double s) { return a *= s; }
  friend constexpr Vec3 operator*(double s, Vec3 a) { return a *= s; }
  friend constexpr Vec3 operator/(Vec3 a, double s) { return a /= s; }
};

[[nodiscard]] constexpr double dot(const Vec3& a, const Vec3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
[[nodiscard]] constexpr Vec3 cross(const Vec3& a, const Vec3& b) {
  return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
[[nodiscard]] inline double norm(const Vec3& a) { return std::sqrt(dot(a, a)); }
[[nodiscard]] inline double distance(const Vec3& a, const Vec3& b) { return norm(a - b); }

/// Returns a unit vector; throws GeometryError for (near) zero vectors.
[[nodiscard]] Vec3 normalized(const Vec3& a);

/// Homogeneous (weighted) point: (w*x, w*y, w*z, w).
struct Vec4 {
  double x{0.0}, y{0.0}, z{0.0}, w{0.0};
  constexpr Vec4() = default;
  constexpr Vec4(double x_, double y_, double z_, double w_) : x(x_), y(y_), z(z_), w(w_) {}
  static constexpr Vec4 from_point(const Vec3& p, double weight) {
    return {p.x * weight, p.y * weight, p.z * weight, weight};
  }
  [[nodiscard]] constexpr Vec3 xyz() const { return {x, y, z}; }
  [[nodiscard]] Vec3 project() const { return Vec3{x / w, y / w, z / w}; }

  constexpr Vec4& operator+=(const Vec4& o) { x += o.x; y += o.y; z += o.z; w += o.w; return *this; }
  friend constexpr Vec4 operator+(Vec4 a, const Vec4& b) { return a += b; }
  friend constexpr Vec4 operator*(double s, const Vec4& a) { return {s * a.x, s * a.y, s * a.z, s * a.w}; }
};

std::ostream& operator<<(std::ostream& os, const Vec3& v);

}  // namespace nf
