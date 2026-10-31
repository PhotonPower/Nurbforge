#pragma once
/// @file join.hpp
/// Joining curves and surfaces end to end.

#include "nurbforge/core/tolerance.hpp"
#include "nurbforge/nurbs/curve.hpp"
#include "nurbforge/nurbs/surface.hpp"

namespace nf {

/// Joins curve b to the end of curve a. The end point of a and the start point of b must agree
/// within `tolerance`; a remaining gap is closed at its midpoint. Degrees are unified by
/// elevation, b's weights are rescaled to match at the joint, and b's parameter interval is
/// shifted to start at a's domain end (keeping its length). The joint knot is then removed as far
/// as this is exact (tol::kLength), so joining the two halves of a split restores the original
/// continuity. Throws GeometryError if the curves do not meet.
[[nodiscard]] NurbsCurve join(const NurbsCurve& a, const NurbsCurve& b,
                              double tolerance = tol::kLength);

/// Joins surface b to surface a across the boundary u = u_max of a / u = u_min of b. The surfaces
/// are made compatible in v (degree elevation, b reparametrised to a's v domain, common knots)
/// and in the u degree. The boundary control rows must agree within `tolerance`, and for
/// rational surfaces their weights must be proportional. b's u interval is shifted to start at
/// a's u_max. The seam knot is removed as far as this is exact. Throws GeometryError otherwise.
[[nodiscard]] NurbsSurface join_u(const NurbsSurface& a, const NurbsSurface& b,
                                  double tolerance = tol::kLength);

/// As join_u, across v = v_max of a / v = v_min of b.
[[nodiscard]] NurbsSurface join_v(const NurbsSurface& a, const NurbsSurface& b,
                                  double tolerance = tol::kLength);

}  // namespace nf
