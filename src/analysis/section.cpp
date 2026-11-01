#include "nurbforge/analysis/section.hpp"

#include <algorithm>
#include <array>
#include <cmath>

#include "nurbforge/core/error.hpp"
#include "nurbforge/geometry/primitives.hpp"
#include "nurbforge/geometry/transform.hpp"

namespace nf {

double Plane::distance(const Vec3& p) const { return dot(p - point, normalized(normal)); }

namespace {

// Orthonormal in-plane axes and a rectangle in the plane covering the projection of `pts`.
struct PlanePatch {
  Vec3 origin, e1, e2;
  double s0{0.0}, s1{1.0}, t0{0.0}, t1{1.0};
  [[nodiscard]] NurbsSurface surface() const {
    return make_rectangle(origin + s0 * e1 + t0 * e2, (s1 - s0) * e1, (t1 - t0) * e2);
  }
  // Maps rectangle parameters (u, v) to in-plane coordinates (s, t).
  [[nodiscard]] std::array<double, 2> st(double u, double v) const {
    return {s0 + u * (s1 - s0), t0 + v * (t1 - t0)};
  }
  [[nodiscard]] Transform uv_to_st() const {
    return Transform({s1 - s0, 0.0, 0.0, s0, 0.0, t1 - t0, 0.0, t0, 0.0, 0.0, 1.0, 0.0});
  }
};

PlanePatch patch_for(const Plane& plane, const std::vector<Vec3>& pts) {
  const Vec3 n = normalized(plane.normal);
  const Vec3 helper = std::abs(n.x) < 0.9 ? Vec3{1, 0, 0} : Vec3{0, 1, 0};
  PlanePatch p{plane.point, normalized(cross(n, helper)), {}};
  p.e2 = cross(n, p.e1);
  double s0 = 1e300, s1 = -1e300, t0 = 1e300, t1 = -1e300;
  for (const auto& q : pts) {
    const double s = dot(q - plane.point, p.e1), t = dot(q - plane.point, p.e2);
    s0 = std::min(s0, s);
    s1 = std::max(s1, s);
    t0 = std::min(t0, t);
    t1 = std::max(t1, t);
  }
  // margin relative to the whole point set (a curve normal to the plane projects to one point)
  Vec3 lo = pts.front(), hi = pts.front();
  for (const auto& q : pts) {
    lo = {std::min(lo.x, q.x), std::min(lo.y, q.y), std::min(lo.z, q.z)};
    hi = {std::max(hi.x, q.x), std::max(hi.y, q.y), std::max(hi.z, q.z)};
  }
  const double margin = 0.1 * std::max({s1 - s0, t1 - t0, distance(lo, hi)}) + 1e-9;
  p.s0 = s0 - margin;
  p.s1 = s1 + margin;
  p.t0 = t0 - margin;
  p.t1 = t1 + margin;
  return p;
}

}  // namespace

SurfaceIntersection planar_section(const NurbsSurface& s, const Plane& plane,
                                   const SurfaceIntersectionOptions& options) {
  const PlanePatch patch = patch_for(plane, s.control_points());
  SurfaceIntersection r = intersect(s, patch.surface(), options);
  const Transform map = patch.uv_to_st();
  for (auto& c : r.curves) {
    c.pcurve2 = c.pcurve2.transformed(map);
    for (auto& p : c.points) {
      const auto st = patch.st(p.u2, p.v2);
      p.u2 = st[0];
      p.v2 = st[1];
    }
  }
  for (auto& p : r.points) {
    const auto st = patch.st(p.u2, p.v2);
    p.u2 = st[0];
    p.v2 = st[1];
  }
  return r;
}

std::vector<SurfaceIntersection> planar_sections(const NurbsSurface& s, const Plane& plane, double spacing,
                                                 std::size_t count, const SurfaceIntersectionOptions& options) {
  if (!(spacing > 0.0)) throw GeometryError("planar_sections: spacing must be positive");
  std::vector<SurfaceIntersection> out;
  const Vec3 n = normalized(plane.normal);
  for (std::size_t k = 0; k < count; ++k)
    out.push_back(planar_section(s, {plane.point + (static_cast<double>(k) * spacing) * n, n}, options));
  return out;
}

std::vector<CurveSurfaceIntersection> intersect(const NurbsCurve& c, const Plane& plane, const Tolerance& tol) {
  const PlanePatch patch = patch_for(plane, c.control_points());
  auto hits = intersect(c, patch.surface(), tol);
  for (auto& h : hits) {
    const auto st = patch.st(h.u, h.v);
    h.u = st[0];
    h.v = st[1];
  }
  return hits;
}

}  // namespace nf
