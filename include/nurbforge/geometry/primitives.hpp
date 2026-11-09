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

/// Full circle through three non-collinear points, starting at p0 and running through p1 and p2
/// in this order (exact, rational). Throws GeometryError for (nearly) collinear points.
[[nodiscard]] NurbsCurve make_circle_3pt(const Vec3& p0, const Vec3& p1, const Vec3& p2);

/// Circular arc from p0 through p1 to p2 (exact, rational; any sweep below 360 degrees).
/// Throws GeometryError for (nearly) collinear points.
[[nodiscard]] NurbsCurve make_arc_3pt(const Vec3& p0, const Vec3& p1, const Vec3& p2);

/// Closed regular polygon with n >= 3 vertices on the circle of the given (circum)radius about
/// `center` in the plane with normal `normal`, counter-clockwise about the normal; the first
/// vertex lies in the direction of x_axis (projected into the plane; zero: any_perpendicular).
[[nodiscard]] NurbsCurve make_polygon(const Vec3& center, const Vec3& normal, double radius,
                                      std::size_t n, const Vec3& x_axis = {0, 0, 0});

/// Closed rectangle outline origin -> origin + edge_u -> + edge_v -> origin + edge_v -> origin
/// (degree-1 curve; the boundary of make_rectangle).
[[nodiscard]] NurbsCurve make_rectangle_curve(const Vec3& origin, const Vec3& edge_u,
                                              const Vec3& edge_v);

/// Helix or spiral about the axis through `axis_point` with direction `axis`, starting at
/// `start`: the radius changes linearly from the start radius (distance of `start` from the
/// axis) to `end_radius` (< 0: the start radius, a cylindrical helix), the height along the axis
/// grows by `pitch` per turn (0: a planar Archimedean spiral; negative: downwards), and the curve
/// winds `turns` > 0 times counter-clockwise about the axis direction (right-handed for
/// pitch > 0), or clockwise if `left_handed`. Conical helices: pitch != 0 and end_radius !=
/// start radius.
///
/// The curve parameter is the number of turns, t in [0, turns]. A helix is not a NURBS; the result
/// is a C2 quintic B-spline (piecewise quintic Hermite interpolation of position, first and second
/// derivative of the exact curve at equally spaced parameters, interior knots of multiplicity 3)
/// whose distance from the exact point at the same parameter is verified below `tolerance` at
/// 24 points per span; the span count is increased (from the observed n^-6 convergence) until
/// it is.
[[nodiscard]] NurbsCurve make_helix(const Vec3& axis_point, const Vec3& axis, const Vec3& start,
                                    double pitch, double turns, double end_radius = -1.0,
                                    bool left_handed = false, double tolerance = 1e-7);

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
