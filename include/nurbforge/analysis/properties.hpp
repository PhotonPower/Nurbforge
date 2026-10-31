#pragma once
/// @file properties.hpp
/// Differential and integral properties: curvature, torsion, arc length, surface area.

#include "nurbforge/core/vec.hpp"
#include "nurbforge/nurbs/curve.hpp"
#include "nurbforge/nurbs/surface.hpp"

namespace nf {

/// Default relative accuracy of the adaptive integrators.
inline constexpr double kIntegrationRelTol = 1e-10;

// ------------------------------------------------------------------ curves

/// Curvature kappa = |C' x C''| / |C'|^3 (>= 0). Throws GeometryError where C' vanishes.
[[nodiscard]] double curvature(const NurbsCurve& c, double u);

/// Torsion tau = (C' x C'') . C''' / |C' x C''|^2. Zero for planar curves; throws GeometryError
/// where the curvature vanishes (torsion undefined).
[[nodiscard]] double torsion(const NurbsCurve& c, double u);

/// Arc length of the curve between the parameters u0 and u1 (u0 > u1 gives a negative value).
/// Adaptive Gauss-Kronrod (7/15) on every knot span, relative accuracy `rel_tol`.
[[nodiscard]] double arc_length(const NurbsCurve& c, double u0, double u1,
                                double rel_tol = kIntegrationRelTol);

/// Arc length of the whole curve.
[[nodiscard]] double arc_length(const NurbsCurve& c, double rel_tol = kIntegrationRelTol);

/// Parameter u with arc_length(c, domain_start, u) = s, for 0 <= s <= total length
/// (safeguarded Newton on the arc length function). Throws GeometryError for s out of range.
[[nodiscard]] double parameter_at_length(const NurbsCurve& c, double s,
                                         double rel_tol = kIntegrationRelTol);

// ------------------------------------------------------------------ surfaces

/// Curvature of a surface at a point, relative to the unit normal n = S_u x S_v / |S_u x S_v|.
struct SurfaceCurvature {
  double gaussian{0.0};  ///< K = k1 * k2
  double mean{0.0};      ///< H = (k1 + k2) / 2 (sign depends on the normal orientation)
  double k1{0.0};        ///< maximal principal curvature
  double k2{0.0};        ///< minimal principal curvature
  Vec3 normal;           ///< unit normal used for the signs
  Vec3 dir1;             ///< unit principal direction of k1 (tangent vector)
  Vec3 dir2;             ///< unit principal direction of k2 (tangent vector)
};

/// Gaussian, mean and principal curvatures and directions from the first and second fundamental
/// forms. At singular points of the parametrisation (poles, collapsed edges) the values are
/// extrapolated (Richardson, from offsets 1e-4 and 2e-4 of the domain towards the interior;
/// typically accurate to ~1e-8 relative for geometrically smooth surfaces). Umbilic points
/// (k1 = k2) return dir1 along S_u and dir2 perpendicular to it.
[[nodiscard]] SurfaceCurvature curvature(const NurbsSurface& s, double u, double v);

/// Area of the surface over [u0, u1] x [v0, v1] (adaptive tensor Gauss-Kronrod on every pair of
/// knot spans, relative accuracy `rel_tol`).
[[nodiscard]] double area(const NurbsSurface& s, double u0, double u1, double v0, double v1,
                          double rel_tol = kIntegrationRelTol);

/// Area of the whole surface.
[[nodiscard]] double area(const NurbsSurface& s, double rel_tol = kIntegrationRelTol);

}  // namespace nf
