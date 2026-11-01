#pragma once
/// @file offset.hpp
/// Offset curves (planar) and offset surfaces, fitted as cubic B-splines to a verified tolerance.
///
/// The exact offset of a NURBS is in general not a NURBS (the unit normal involves a square
/// root), so the offset is approximated by a non-rational cubic spline. The parametrisation of
/// the input is kept: offset(t) approximates c(t) + d n(t) (resp. S(u,v) + d N(u,v)), and the
/// distance between the two is checked against `tolerance` at dense points of every fitted span;
/// spans are refined adaptively until the check passes.
///
/// Only the regular case is handled: an offset towards the centre of curvature by at least the
/// radius of curvature (cusps / degenerate normals) is rejected with GeometryError. Global
/// self-intersections of the offset (e.g. a bottleneck narrower than 2|d|) are not detected and
/// not removed.

#include "nurbforge/core/vec.hpp"
#include "nurbforge/nurbs/curve.hpp"
#include "nurbforge/nurbs/surface.hpp"

namespace nf {

/// Planar offset c(t) + d * normalize(plane_normal x c'(t)) of a curve lying in a plane with the
/// given normal (d > 0: to the left of the direction of travel, seen from the tip of
/// plane_normal; e.g. inwards for a counter-clockwise circle about +z).
///
/// The result is a non-rational cubic B-spline on the domain of c with offset(t) belonging to
/// c(t). Every Bézier segment of c is fitted separately (rational curves such as circles are only
/// C1 at their knots) to tolerance/2 and the pieces are joined with nf::join, so the result has
/// (at most) C0 knots at the knots of c. Throws GeometryError if the control points do not lie in
/// a plane normal to plane_normal (within `tolerance`), if c' vanishes, if c has a tangent
/// discontinuity at a knot (the offset would have a gap larger than tolerance/2), if
/// d * signed curvature >= 1 somewhere (cusp: |d| reaches the radius of curvature on the concave
/// side), or if the tolerance cannot be reached.
[[nodiscard]] NurbsCurve offset_curve(const NurbsCurve& c, double distance,
                                      const Vec3& plane_normal, double tolerance = 1e-7);

/// Offset surface S(u,v) + d * N(u,v) with the unit normal N = normal(u, v) (limit normals at
/// poles, so offsets of revolved surfaces with poles are well defined).
///
/// The result is a non-rational tensor-product spline of degree 5 x 5 (quintic: the error
/// decreases like h^6, so far fewer control points are needed than with a bicubic at tight
/// tolerances) on the domain of s with offset(u,v) belonging to S(u,v). The knot lines of s are
/// kept as break lines (with a multiplicity matching the continuity of s, up to C0). Knots are
/// first inserted adaptively from fits along iso-lines (to tolerance/2 per direction), then the
/// tensor-product fit is verified at 9 x 9 check points of every knot rectangle (failing
/// rectangles are refined further). Exact
/// references: plane -> plane, sphere -> concentric sphere R + d, cylinder -> coaxial cylinder,
/// torus -> torus with tube radius r + d (signs follow the normal orientation).
/// Throws GeometryError if d * k >= 1 for a principal curvature k (relative to N) on a sample
/// grid (9 x 9 points per knot rectangle), if the normal is discontinuous across a knot line by
/// more than the tolerance allows, or if the tolerance cannot be reached. Local self-intersections
/// are not removed.
[[nodiscard]] NurbsSurface offset_surface(const NurbsSurface& s, double distance,
                                          double tolerance = 1e-7);

}  // namespace nf
