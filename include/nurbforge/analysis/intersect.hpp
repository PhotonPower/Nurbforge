#pragma once
/// @file intersect.hpp
/// Intersections of curves with curves and with surfaces.

#include <vector>

#include "nurbforge/core/tolerance.hpp"
#include "nurbforge/core/vec.hpp"
#include "nurbforge/nurbs/curve.hpp"
#include "nurbforge/nurbs/surface.hpp"

namespace nf {

struct CurveCurveIntersection {
  double s{0.0};         ///< parameter on the first curve
  double t{0.0};         ///< parameter on the second curve
  Vec3 point;            ///< midpoint of a(s) and b(t)
  double distance{0.0};  ///< |a(s) - b(t)| (<= tol.length)
};

struct CurveSurfaceIntersection {
  double t{0.0};         ///< parameter on the curve
  double u{0.0}, v{0.0}; ///< parameters on the surface
  Vec3 point;            ///< midpoint of c(t) and S(u, v)
  double distance{0.0};  ///< |c(t) - S(u, v)| (<= tol.length)
};

/// All intersection points of two curves (in 3D: points where they come within tol.length).
/// Subdivision of Bézier pieces with bounding-box rejection until both pieces are monotone,
/// then box-constrained Gauss-Newton on a(s) - b(t) = 0; duplicates at piece boundaries are
/// merged. Sorted by s. Tangential intersections are found as long as the curves come within
/// tol.length; coincident (overlapping) curve portions are not detected as such and yield
/// sample points only.
[[nodiscard]] std::vector<CurveCurveIntersection> intersect(const NurbsCurve& a, const NurbsCurve& b,
                                                            const Tolerance& tol = {});

/// All intersection points of a curve and a surface, sorted by t (same method, Gauss-Newton on
/// c(t) - S(u, v) = 0). Curve portions lying in the surface are not detected as such.
[[nodiscard]] std::vector<CurveSurfaceIntersection> intersect(const NurbsCurve& c, const NurbsSurface& s,
                                                              const Tolerance& tol = {});

}  // namespace nf
