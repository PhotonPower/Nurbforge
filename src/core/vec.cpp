#include "nurbforge/core/vec.hpp"
#include "nurbforge/core/error.hpp"

namespace nf {

Vec3 normalized(const Vec3& a) {
  const double n = norm(a);
  if (n < kTolerance) throw GeometryError("cannot normalize a zero-length vector");
  return a / n;
}

std::ostream& operator<<(std::ostream& os, const Vec3& v) {
  return os << '(' << v.x << ", " << v.y << ", " << v.z << ')';
}

}  // namespace nf
