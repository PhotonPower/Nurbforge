#include "nurbforge/geometry/transform.hpp"

#include <cmath>

namespace nf {

Transform::Transform() : m_{1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0} {}

Transform Transform::translation(const Vec3& d) {
  return Transform({1, 0, 0, d.x, 0, 1, 0, d.y, 0, 0, 1, d.z});
}

Transform Transform::rotation(const Vec3& origin, const Vec3& axis, double angle) {
  const Vec3 k = normalized(axis);
  const double c = std::cos(angle), s = std::sin(angle), t = 1.0 - c;
  // Rodrigues rotation matrix
  const std::array<double, 9> R{
      t * k.x * k.x + c,       t * k.x * k.y - s * k.z, t * k.x * k.z + s * k.y,
      t * k.x * k.y + s * k.z, t * k.y * k.y + c,       t * k.y * k.z - s * k.x,
      t * k.x * k.z - s * k.y, t * k.y * k.z + s * k.x, t * k.z * k.z + c};
  const Transform rot({R[0], R[1], R[2], 0, R[3], R[4], R[5], 0, R[6], R[7], R[8], 0});
  return translation(origin) * rot * translation(-origin);
}

Transform Transform::scaling(const Vec3& center, double f) { return scaling(center, {f, f, f}); }

Transform Transform::scaling(const Vec3& center, const Vec3& f) {
  const Transform s({f.x, 0, 0, 0, 0, f.y, 0, 0, 0, 0, f.z, 0});
  return translation(center) * s * translation(-center);
}

Transform Transform::mirror(const Vec3& point, const Vec3& normal) {
  const Vec3 n = normalized(normal);
  // Householder: I - 2 n n^T
  const Transform h({1 - 2 * n.x * n.x, -2 * n.x * n.y, -2 * n.x * n.z, 0,
                     -2 * n.y * n.x, 1 - 2 * n.y * n.y, -2 * n.y * n.z, 0,
                     -2 * n.z * n.x, -2 * n.z * n.y, 1 - 2 * n.z * n.z, 0});
  return translation(point) * h * translation(-point);
}

Vec3 Transform::apply(const Vec3& p) const {
  return {m_[0] * p.x + m_[1] * p.y + m_[2] * p.z + m_[3],
          m_[4] * p.x + m_[5] * p.y + m_[6] * p.z + m_[7],
          m_[8] * p.x + m_[9] * p.y + m_[10] * p.z + m_[11]};
}

Vec3 Transform::apply_vector(const Vec3& v) const {
  return {m_[0] * v.x + m_[1] * v.y + m_[2] * v.z,
          m_[4] * v.x + m_[5] * v.y + m_[6] * v.z,
          m_[8] * v.x + m_[9] * v.y + m_[10] * v.z};
}

Transform operator*(const Transform& a, const Transform& b) {
  const auto& A = a.m_;
  const auto& B = b.m_;
  std::array<double, 12> C{};
  for (int r = 0; r < 3; ++r) {
    for (int c = 0; c < 4; ++c) {
      double v = 0.0;
      for (int k = 0; k < 3; ++k) v += A[static_cast<std::size_t>(4 * r + k)] * B[static_cast<std::size_t>(4 * k + c)];
      if (c == 3) v += A[static_cast<std::size_t>(4 * r + 3)];
      C[static_cast<std::size_t>(4 * r + c)] = v;
    }
  }
  return Transform(C);
}

}  // namespace nf
