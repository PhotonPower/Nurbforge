#pragma once
/// @file surface_intersect.hpp
/// Surface/surface intersection: seeds by subdivision, curves by marching.

#include <vector>

#include "nurbforge/core/tolerance.hpp"
#include "nurbforge/core/vec.hpp"
#include "nurbforge/nurbs/curve.hpp"
#include "nurbforge/nurbs/surface.hpp"

namespace nf {

/// A point on both surfaces with its parameters.
struct SurfaceSurfacePoint {
  Vec3 point;
  double u1{0.0}, v1{0.0};  ///< parameters on the first surface
  double u2{0.0}, v2{0.0};  ///< parameters on the second surface
};

/// One intersection curve.
struct IntersectionCurve {
  std::vector<SurfaceSurfacePoint> points;  ///< traced points in order (first == last if closed)
  bool closed{false};
  NurbsCurve curve;    ///< cubic interpolation of the points in 3D (chord-length parameter)
  NurbsCurve pcurve1;  ///< (u1, v1, 0) in the parameter domain of the first surface, same parameter
  NurbsCurve pcurve2;  ///< (u2, v2, 0) in the parameter domain of the second surface
};

struct SurfaceIntersection {
  std::vector<IntersectionCurve> curves;
  std::vector<SurfaceSurfacePoint> points;  ///< isolated (tangential) contact points
};

struct SurfaceIntersectionOptions {
  Tolerance tol{};
  /// Maximal deviation of the traced polyline from the true curve (model units);
  /// 0: 1e-4 of the model size. Controls the marching step.
  double chord_tolerance{0.0};
  /// Accuracy of the fitted curves (model units): |S1(pcurve1(t)) - curve(t)| and
  /// |S2(pcurve2(t)) - curve(t)| are verified between all fit points and exact points are added
  /// where needed (e.g. near poles, where the p-curves vary quickly); 0: 1e-7 of the model size.
  double fit_tolerance{0.0};
};

/// Intersection of two surfaces. Seeds are found by subdividing pairs of Bézier patches; a pair
/// whose normal cones are separated contains no closed loop, so all its intersection curves
/// cross the pair's boundary and are found by intersecting the boundary curves with the other
/// patch (curve/surface intersection). Curves are traced by predictor-corrector marching along
/// n1 x n2 with step control by the chord tolerance; they end on the boundary of either
/// parameter domain (solved exactly there) or close into loops. Curves are split where they
/// cross the domain boundary (e.g. the seam of a closed surface). Points where the surfaces
/// only touch are reported in `points`; surfaces touching along a curve (tangential
/// intersection curves) are not traced.
[[nodiscard]] SurfaceIntersection intersect(const NurbsSurface& a, const NurbsSurface& b,
                                            const SurfaceIntersectionOptions& options = {});

}  // namespace nf
