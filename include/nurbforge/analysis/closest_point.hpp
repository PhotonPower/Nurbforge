#pragma once
/// @file closest_point.hpp
/// Closest point (point projection / point inversion) on curves and surfaces.

#include "nurbforge/core/tolerance.hpp"
#include "nurbforge/core/vec.hpp"
#include "nurbforge/nurbs/curve.hpp"
#include "nurbforge/nurbs/surface.hpp"

namespace nf {

/// Foot point on a curve.
struct CurvePoint {
  double u{0.0};        ///< parameter of the foot point
  Vec3 point;           ///< C(u)
  double distance{0.0}; ///< |C(u) - P|
};

/// Foot point on a surface.
struct SurfacePoint {
  double u{0.0}, v{0.0};
  Vec3 point;           ///< S(u, v)
  double distance{0.0}; ///< |S(u, v) - P|
};

/// Global closest point of P on the curve. Branch and bound over Bézier segments (the distance
/// to the bounding box of a segment's control points is a safe lower bound, so the global
/// minimum is never discarded), followed by safeguarded Newton iteration (P&T Eq. 6.3) in every
/// remaining monotone segment. Uses tol.length (point coincidence, step size) and tol.angular
/// (orthogonality of C' and C - P). If several points are equally close (e.g. the centre of a
/// circle), one of them is returned.
[[nodiscard]] CurvePoint closest_point(const NurbsCurve& c, const Vec3& P, const Tolerance& tol = {});

/// Local Newton projection starting at u0 (domain-clamped). Fast, but only finds the local
/// minimum near u0; intended for marching and refinement of a known approximation.
[[nodiscard]] CurvePoint closest_point_local(const NurbsCurve& c, const Vec3& P, double u0,
                                             const Tolerance& tol = {});

/// Global closest point of P on the surface (branch and bound over Bézier patches, then
/// safeguarded 2D Newton, P&T Eq. 6.6). Robust at poles and on the domain boundary.
[[nodiscard]] SurfacePoint closest_point(const NurbsSurface& s, const Vec3& P,
                                         const Tolerance& tol = {});

/// Local 2D Newton projection starting at (u0, v0).
[[nodiscard]] SurfacePoint closest_point_local(const NurbsSurface& s, const Vec3& P, double u0,
                                               double v0, const Tolerance& tol = {});

}  // namespace nf
