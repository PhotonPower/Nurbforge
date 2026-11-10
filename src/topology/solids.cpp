#include "nurbforge/topology/solids.hpp"

#include <cmath>
#include <numbers>
#include <string>
#include <utility>
#include <vector>

#include "nurbforge/core/error.hpp"
#include "nurbforge/geometry/primitives.hpp"
#include "nurbforge/geometry/transform.hpp"
#include "nurbforge/modeling/construct.hpp"

namespace nf {

namespace {

constexpr double kPi = std::numbers::pi;

/// Map from the local frame (z along `axis`) to the world, origin at `origin`.
Transform frame(const Vec3& origin, const Vec3& axis) {
  const double l = norm(axis);
  if (!(l > 0.0)) throw GeometryError("solid primitive: axis must not be zero");
  const Vec3 z = (1.0 / l) * axis;
  const Vec3 x =
      normalized(std::abs(z.x) < 0.9 ? cross(z, Vec3{1, 0, 0}) : cross(z, Vec3{0, 1, 0}));
  const Vec3 y = cross(z, x);
  return Transform({x.x, y.x, z.x, origin.x, x.y, y.y, z.y, origin.y, x.z, y.z, z.z, origin.z});
}

/// Planar disc of radius r at height z (local frame): square trimmed by an exact circle.
TrimmedSurface disc(double r, double z) {
  return TrimmedSurface(make_rectangle({-r, -r, z}, {2 * r, 0, 0}, {0, 2 * r, 0}),
                        circle_loop(0.5, 0.5, 0.5));
}

/// Sews faces built in the local frame after mapping them by t, outward oriented.
Shell assemble(std::vector<TrimmedSurface> local, const Transform& t) {
  std::vector<TrimmedSurface> faces;
  faces.reserve(local.size());
  for (const auto& f : local) faces.emplace_back(f.surface().transformed(t), f.outer(), f.holes());
  return Shell::sew(std::move(faces)).oriented();
}

void require_positive(double v, const char* what) {
  if (!(v > 0.0) || !std::isfinite(v))
    throw GeometryError(std::string("solid primitive: ") + what + " must be > 0");
}

}  // namespace

Shell make_box(const Vec3& corner, const Vec3& size) {
  require_positive(size.x, "box size x");
  require_positive(size.y, "box size y");
  require_positive(size.z, "box size z");
  const Vec3 c = corner, ex{size.x, 0, 0}, ey{0, size.y, 0}, ez{0, 0, size.z};
  return Shell::sew({TrimmedSurface(make_rectangle(c, ex, ey)),
                     TrimmedSurface(make_rectangle(c + ez, ex, ey)),
                     TrimmedSurface(make_rectangle(c, ex, ez)),
                     TrimmedSurface(make_rectangle(c + ey, ex, ez)),
                     TrimmedSurface(make_rectangle(c, ey, ez)),
                     TrimmedSurface(make_rectangle(c + ex, ey, ez))})
      .oriented();
}

Shell make_sphere(const Vec3& center, double radius) {
  require_positive(radius, "sphere radius");
  const auto meridian = make_arc(center, {0, 0, -1}, {1, 0, 0}, radius, 0.0, kPi);
  return Shell::sew({TrimmedSurface(revolve(meridian, center, {0, 0, 1}, 2 * kPi))}).oriented();
}

Shell make_cylinder(const Vec3& base, const Vec3& axis, double radius, double height) {
  require_positive(radius, "cylinder radius");
  require_positive(height, "cylinder height");
  const auto side = extrude(make_circle({0, 0, 0}, {0, 0, 1}, radius), {0, 0, height});
  return assemble({TrimmedSurface(side), disc(radius, 0.0), disc(radius, height)},
                  frame(base, axis));
}

Shell make_cone(const Vec3& base, const Vec3& axis, double r0, double r1, double height) {
  require_positive(height, "cone height");
  if (!(r0 >= 0.0) || !(r1 >= 0.0) || !(r0 > 0.0 || r1 > 0.0))
    throw GeometryError("solid primitive: cone radii must be >= 0 and not both 0");
  const auto side = revolve(make_line({r0, 0, 0}, {r1, 0, height}), {0, 0, 0}, {0, 0, 1}, 2 * kPi);
  std::vector<TrimmedSurface> faces{TrimmedSurface(side)};
  if (r0 > 0) faces.push_back(disc(r0, 0.0));
  if (r1 > 0) faces.push_back(disc(r1, height));
  return assemble(std::move(faces), frame(base, axis));
}

Shell make_torus(const Vec3& center, const Vec3& axis, double major_radius, double minor_radius) {
  require_positive(minor_radius, "torus minor radius");
  if (!(major_radius > minor_radius))
    throw GeometryError("solid primitive: torus major radius must exceed the minor radius");
  const auto tube = revolve(make_circle({major_radius, 0, 0}, {0, 1, 0}, minor_radius), {0, 0, 0},
                            {0, 0, 1}, 2 * kPi);
  return assemble({TrimmedSurface(tube)}, frame(center, axis));
}

}  // namespace nf
