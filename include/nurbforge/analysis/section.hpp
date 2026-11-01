#pragma once
/// @file section.hpp
/// Planar sections of curves and surfaces.

#include <cstddef>
#include <vector>

#include "nurbforge/analysis/intersect.hpp"
#include "nurbforge/analysis/surface_intersect.hpp"

namespace nf {

/// Infinite plane through `point` with normal `normal` (normalised internally).
struct Plane {
  Vec3 point;
  Vec3 normal{0, 0, 1};
  /// Signed distance of p from the plane.
  [[nodiscard]] double distance(const Vec3& p) const;
};

/// Section of a surface with a plane: the intersection curves (3D curve and p-curve on the
/// surface; pcurve2 is the curve in the plane's own (s, t) coordinates along two orthonormal
/// in-plane axes, origin at plane.point) and touching points. The plane is represented by a
/// rectangle covering the surface's control-point bounding box, so all section curves are
/// found; curves are split where they cross the surface's parameter boundary (seams).
[[nodiscard]] SurfaceIntersection planar_section(const NurbsSurface& s, const Plane& plane,
                                                 const SurfaceIntersectionOptions& options = {});

/// `count` parallel sections, the k-th at plane.point + k * spacing * normal (k = 0..count-1).
[[nodiscard]] std::vector<SurfaceIntersection> planar_sections(const NurbsSurface& s, const Plane& plane,
                                                               double spacing, std::size_t count,
                                                               const SurfaceIntersectionOptions& options = {});

/// Intersection points of a curve with a plane (curve parameter t, point; u, v are the in-plane
/// coordinates as for planar_section).
[[nodiscard]] std::vector<CurveSurfaceIntersection> intersect(const NurbsCurve& c, const Plane& plane,
                                                              const Tolerance& tol = {});

}  // namespace nf
