#pragma once
/// @file blend.hpp
/// Blend curves (smooth transitions between the ends of two curves), tween curves, and lines
/// tangent to curves.

#include <cstddef>
#include <vector>

#include "nurbforge/modeling/extend.hpp"
#include "nurbforge/nurbs/curve.hpp"

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

}  // namespace nf
