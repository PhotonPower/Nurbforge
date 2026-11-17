#pragma once
/// @file blend.hpp
/// Blend curves: smooth transitions between the ends of two curves.

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

}  // namespace nf
