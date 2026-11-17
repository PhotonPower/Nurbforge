#pragma once
/// @file project.hpp
/// Projection of curves onto surfaces.

#include "nurbforge/nurbs/curve.hpp"
#include "nurbforge/nurbs/surface.hpp"

namespace nf {

/// A curve on a surface: the 3D curve and its p-curve (u, v, 0) in the surface's parameter
/// domain, both with the parameter of the projected curve.
struct ProjectedCurve {
  NurbsCurve curve;   ///< 3D curve on the surface
  NurbsCurve pcurve;  ///< (u, v, 0) on the surface; S(pcurve(t)) = curve(t) within the tolerance
};

/// Projects a curve onto a surface: along the surface normal (closest points; `direction` zero)
/// or along `direction` (every point moved parallel to it until it meets the surface). The
/// projected points are interpolated (cubic, with the parameters of the input curve) and points
/// are added until the 3D curve and the image of the p-curve are within `tolerance` of the
/// exact projections at the midpoints between them. Throws GeometryError if a point does not
/// project (a ray misses the surface), or if the projection crosses the seam of a closed
/// surface (split the curve there).
[[nodiscard]] ProjectedCurve project_curve(const NurbsCurve& curve, const NurbsSurface& surface,
                                           const Vec3& direction = {0, 0, 0},
                                           double tolerance = 1e-6);

}  // namespace nf
