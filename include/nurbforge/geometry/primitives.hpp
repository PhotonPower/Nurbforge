#pragma once
/// @file primitives.hpp
/// Exact NURBS representations of elementary geometry.

#include "nurbforge/core/vec.hpp"
#include "nurbforge/nurbs/curve.hpp"
#include "nurbforge/nurbs/surface.hpp"

namespace nf {

/// Straight line segment p0 -> p1 (degree 1).
[[nodiscard]] NurbsCurve make_line(const Vec3& p0, const Vec3& p1);

/// Open polyline through the given points (degree 1).
[[nodiscard]] NurbsCurve make_polyline(const std::vector<Vec3>& points);

/// Circular arc (A7.1). x_axis and y_axis span the arc plane and must be orthogonal.
/// Angles in radians, measured from x_axis towards y_axis.
[[nodiscard]] NurbsCurve make_arc(const Vec3& center, const Vec3& x_axis, const Vec3& y_axis,
                                  double radius, double start_angle, double end_angle);

/// Full circle in the plane through `center` with the given normal (exact, 9 control points).
[[nodiscard]] NurbsCurve make_circle(const Vec3& center, const Vec3& normal, double radius);

/// Bilinear patch from 4 corners: S(0,0)=p00, S(1,0)=p10, S(0,1)=p01, S(1,1)=p11.
[[nodiscard]] NurbsSurface make_bilinear(const Vec3& p00, const Vec3& p10, const Vec3& p01,
                                         const Vec3& p11);

/// Planar rectangle with corner `origin`, edges `edge_u` and `edge_v`.
[[nodiscard]] NurbsSurface make_rectangle(const Vec3& origin, const Vec3& edge_u, const Vec3& edge_v);

/// Any unit vector orthogonal to n (deterministic choice).
[[nodiscard]] Vec3 any_perpendicular(const Vec3& n);

namespace detail {
/// Number of arcs (<= 90 deg each) used for an angular sweep, and the matching knot vector.
std::size_t arc_segment_count(double sweep);
KnotVector arc_knots(std::size_t narcs);
/// Intersection of two coplanar lines P0 + a*T0 and P2 + b*T2.
Vec3 intersect_lines(const Vec3& p0, const Vec3& t0, const Vec3& p2, const Vec3& t2);
}  // namespace detail

}  // namespace nf
