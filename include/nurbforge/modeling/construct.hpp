#pragma once
/// @file construct.hpp
/// Surface construction from curves (extrusion, revolution, ruled surfaces).

#include "nurbforge/nurbs/curve.hpp"
#include "nurbforge/nurbs/surface.hpp"

namespace nf {

/// General cylinder: sweeps `profile` along the straight vector `direction`.
/// Result: u follows the profile, v the extrusion (degree 1).
[[nodiscard]] NurbsSurface extrude(const NurbsCurve& profile, const Vec3& direction);

/// Surface of revolution (A8.1): rotates `profile` by `angle` (radians, 0 < angle <= 2*pi)
/// about the axis through `axis_point` with direction `axis_dir`.
/// Result: u runs around the axis (degree 2, rational), v follows the profile.
[[nodiscard]] NurbsSurface revolve(const NurbsCurve& profile, const Vec3& axis_point,
                                   const Vec3& axis_dir, double angle);

/// Ruled (linearly lofted) surface between two curves.
/// v = 0 lies on c0, v = 1 on c1. Curves are made compatible (degree and knots)
/// automatically, see make_compatible.
[[nodiscard]] NurbsSurface ruled(const NurbsCurve& c0, const NurbsCurve& c1);

/// Makes two curves compatible: clamped knot vectors, common domain [0, 1], equal degree
/// (by elevating the lower one) and a common knot vector (by knot insertion). Shapes are
/// unchanged.
void make_compatible(NurbsCurve& a, NurbsCurve& b);

/// Linear reparametrisation of a curve to the domain [a, b].
[[nodiscard]] NurbsCurve reparametrize(const NurbsCurve& c, double a, double b);

}  // namespace nf
