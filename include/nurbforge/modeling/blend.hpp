#pragma once
/// @file blend.hpp
/// Blend curves (smooth transitions between the ends of two curves), tween curves, lines
/// tangent to curves, and blend / match surfaces (smooth transitions between surface sides).

#include <cstddef>
#include <vector>

#include "nurbforge/modeling/extend.hpp"
#include "nurbforge/nurbs/curve.hpp"
#include "nurbforge/nurbs/surface.hpp"

namespace nf {

/// A Bézier curve from the end `end0` of curve a to the end `end1` of curve b that continues
/// both curves smoothly: G1 (continuity 1, cubic: same tangent directions) or G2 (continuity 2,
/// quintic: same tangent directions and curvature vectors). The blend leaves a away from it
/// (along the tangent of a at End, against it at Start) and enters b towards its interior.
/// The speed |C'| at both ends is `tension` times the chord between the ends; larger values
/// give fuller blends. Domain [0, 1]. Throws GeometryError if the ends coincide or a tangent
/// vanishes.
[[nodiscard]] NurbsCurve blend_curves(const NurbsCurve& a, CurveEnd end0, const NurbsCurve& b,
                                      CurveEnd end1, int continuity = 2, double tension = 1.0);

/// `count` intermediate curves between a (s = 0) and b (s = 1) at s = k / (count + 1),
/// k = 1..count. The curves are made compatible (make_compatible: common degree, knots and
/// domain [0, 1]) and the homogeneous control points (w P, w) are interpolated linearly, so
/// tween k is exactly (1 - s) a + s b in homogeneous space (for non-rational curves: the
/// point-wise affine combination at equal parameters). Both curves should run in the same
/// direction (reverse one otherwise). Throws GeometryError if count is 0.
[[nodiscard]] std::vector<NurbsCurve> tween_curves(const NurbsCurve& a, const NurbsCurve& b,
                                                   std::size_t count);

/// A line segment tangent to curves: from point_a to point_b, touching curve a at param_a (or
/// starting at the given point, param_a NaN) and curve b at param_b.
struct TangentLine {
  NurbsCurve line;  ///< degree 1, from point_a to point_b
  Vec3 point_a, point_b;
  double param_a{0.0}, param_b{0.0};
};

/// Lines through `point` tangent to a planar curve (all solutions, sorted by the curve
/// parameter): the tangent points t with n . ((c(t) - p) x c'(t)) = 0 for the plane normal n,
/// found by sampling and refined by Newton. Each line runs from the point to the tangent point
/// (param_a is NaN). Points on the curve itself are not tangent points. Throws GeometryError if
/// the point or the curve is not in a plane normal to plane_normal (within 1e-9 of the curve
/// size).
[[nodiscard]] std::vector<TangentLine> tangent_lines(const Vec3& point, const NurbsCurve& curve,
                                                     const Vec3& plane_normal = {0, 0, 1});

/// Common tangent lines of two coplanar curves (all solutions found on a sampling grid, refined
/// by Newton): parameters (s, t) where a'(s), b'(t) and b(t) - a(s) are parallel. Each line runs
/// from a(s) to b(t); intersection points of the curves are excluded. Sorted by s, then t.
/// Throws GeometryError if the curves are not in a common plane normal to plane_normal.
[[nodiscard]] std::vector<TangentLine> tangent_lines(const NurbsCurve& a, const NurbsCurve& b,
                                                     const Vec3& plane_normal = {0, 0, 1});

/// Blend surface between the side `side0` of s0 and the side `side1` of s1 that continues both
/// surfaces smoothly: G0 (continuity 0: ruled), G1 (1: same tangent planes along both edges,
/// cubic across) or G2 (2: also the same normal curvature, quintic across). u runs along the
/// edge of s0 (in its parameter direction), v across from s0 (v = 0) to s1 (v = 1); the edge of
/// s1 is mapped linearly onto the edge of s0, in the orientation that pairs the nearer points
/// (for closed edges the seams should be aligned). Across the edges the blend is a
/// reparametrisation of the surfaces up to the requested order: B_v = k D1 and B_vv = k^2 D2
/// with the outward cross-boundary derivatives D1, D2 of the surface and a constant k chosen so
/// that the mean speed |B_v| is `tension` times the mean distance between the edges.
/// The result is exact when the weights of each surface do not vary across its side in the
/// rows involved (non-rational surfaces, extrusions, cylinders along the axis); otherwise the
/// edge and its derivatives are interpolated to `tolerance` (relative to the edge length for the
/// derivatives). Throws GeometryError if the edges touch or a cross derivative vanishes.
[[nodiscard]] NurbsSurface blend_surfaces(const NurbsSurface& s0, SurfaceSide side0,
                                          const NurbsSurface& s1, SurfaceSide side1,
                                          int continuity = 2, double tension = 1.0,
                                          double tolerance = 1e-7);

/// Changes `surface` near its side `side` so that it joins `target` along the side
/// `target_side` with continuity G0 (the edges coincide), G1 (also the tangent planes) or G2
/// (also the normal curvature). The edge of the target is mapped linearly onto the edge of the
/// surface (orientation as in blend_surfaces); the first continuity + 1 control rows at the side
/// are replaced so that the surface continues the target as a reparametrisation (cross speed
/// kept on average). The surface is refined first (degree and knots along the edge, knots
/// across) so that the opposite side and its cross derivatives up to the same order are kept.
/// Exactness as for blend_surfaces. Returns the modified surface with the same parameter
/// directions and domain.
[[nodiscard]] NurbsSurface match_surface(const NurbsSurface& surface, SurfaceSide side,
                                         const NurbsSurface& target, SurfaceSide target_side,
                                         int continuity = 2, double tolerance = 1e-7);

}  // namespace nf
