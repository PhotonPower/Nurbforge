#pragma once
/// @file project.hpp
/// Projection of curves onto surfaces and planes; curves on surfaces (interpolation through
/// points on a surface, offsets within a surface).

#include <cstddef>
#include <vector>

#include "nurbforge/geometry/plane.hpp"
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

/// Parallel projection onto a plane: every point moves along `direction` (zero: along the plane
/// normal, i.e. orthogonal projection) until it meets the plane. The map is affine, so the
/// projection of the control points (weights kept) is exact. Throws GeometryError if the
/// direction is parallel to the plane.
[[nodiscard]] NurbsCurve project_to_plane(const NurbsCurve& curve, const Plane& plane,
                                          const Vec3& direction = {0, 0, 0});
/// Parallel projection of a surface onto a plane (exact, see the curve overload).
[[nodiscard]] NurbsSurface project_to_plane(const NurbsSurface& surface, const Plane& plane,
                                            const Vec3& direction = {0, 0, 0});

/// A curve on a surface through points: the points are projected onto the surface (closest
/// points, continued from the previous point so that the curve stays on one branch), or taken
/// as (u, v, 0) parameter points if `points_are_uv`; the p-curve interpolates them in the
/// parameter domain (degree `degree`, parameters by chord length of the 3D points), and the 3D
/// curve follows S(pcurve(t)) within `tolerance`. Both share the parameter t in [0, 1]. Throws
/// GeometryError for fewer than 2 points, or if consecutive points lie on opposite sides of the
/// seam of a closed surface.
[[nodiscard]] ProjectedCurve interpolate_on_surface(const NurbsSurface& surface,
                                                    const std::vector<Vec3>& points,
                                                    std::size_t degree = 3,
                                                    bool points_are_uv = false,
                                                    double tolerance = 1e-6);

/// Offset of a curve within a surface (geodesic offset): from every point S(pcurve(t)) a
/// geodesic of length |distance| is traced perpendicular to the curve, to the left of its
/// direction seen from the surface normal S_u x S_v for distance > 0 (to the right for
/// distance < 0). The geodesics are integrated with parallel transport of the direction
/// (midpoint rule, steps halved until converged); the end points are interpolated like
/// project_curve until the offset curve and its p-curve are within `tolerance` of the exact
/// offset points. The result shares the parameter of `pcurve`. Throws GeometryError if a
/// geodesic leaves the surface's parameter domain or the surface is degenerate on the way.
[[nodiscard]] ProjectedCurve offset_on_surface(const NurbsSurface& surface,
                                               const NurbsCurve& pcurve, double distance,
                                               double tolerance = 1e-6);

}  // namespace nf
