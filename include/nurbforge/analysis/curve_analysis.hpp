#pragma once
/// @file curve_analysis.hpp
/// Analysis of curves: osculating circle, curvature graph, direction and orientation,
/// continuity between curves and at knots, duplicate curves.

#include <cstddef>
#include <vector>

#include "nurbforge/core/vec.hpp"
#include "nurbforge/modeling/extend.hpp"
#include "nurbforge/nurbs/curve.hpp"

namespace nf {

// ------------------------------------------------------------------ curvature

/// Curvature of a curve at a parameter, with the osculating circle.
struct CurveCurvature {
  Vec3 point;             ///< C(u)
  Vec3 tangent;           ///< unit tangent
  Vec3 normal;            ///< unit principal normal, towards the centre of curvature (zero where
                          ///< the curvature vanishes)
  Vec3 vector;            ///< curvature vector: curvature * normal
  double curvature{0.0};  ///< kappa >= 0
  double radius{0.0};     ///< radius of curvature 1 / kappa (infinity where kappa = 0)
  Vec3 center;            ///< centre of the osculating circle (the point itself where kappa = 0)
};

/// Curvature, radius of curvature and osculating circle at u. Throws GeometryError where C'
/// vanishes.
[[nodiscard]] CurveCurvature curve_curvature(const NurbsCurve& c, double u);

/// Curvature graph (comb) of a curve: at every sample a hair from the curve point to
/// point - scale * curvature vector, i.e. on the convex side, with a length proportional to the
/// curvature.
struct CurvatureGraph {
  std::vector<double> parameters;
  std::vector<Vec3> points;        ///< curve points (feet of the hairs)
  std::vector<Vec3> tips;          ///< ends of the hairs; in order they form the graph line
  std::vector<double> curvatures;  ///< kappa at the samples
  double max_curvature{0.0};
  double max_parameter{0.0};  ///< parameter of the sample with the largest curvature
};

/// Curvature graph with `samples_per_span` + 1 samples on every knot span. The ends of a span
/// are evaluated from within the span, so a knot where the curvature jumps (a line running
/// into an arc: G1 only) appears twice with both values and the graph shows the step.
[[nodiscard]] CurvatureGraph curvature_graph(const NurbsCurve& c, double scale = 1.0,
                                             std::size_t samples_per_span = 16);

// ------------------------------------------------------------------ direction

/// True if the end points of the curve agree within the tolerance.
[[nodiscard]] bool is_closed(const NurbsCurve& c, double tolerance = 1e-9);

/// True if the curve lies in a plane within the tolerance (its control points do).
[[nodiscard]] bool is_planar(const NurbsCurve& c, double tolerance = 1e-9);

/// Direction of a curve: where it starts and ends, and the sense of rotation of a closed curve.
struct CurveDirection {
  Vec3 start, end;                  ///< end points
  Vec3 start_tangent, end_tangent;  ///< unit tangents in the direction of travel
  bool closed{false};
  bool planar{false};
  /// Normal of the plane the sense of rotation refers to: the given one, or the normal of the
  /// curve's plane with its largest component positive (zero for a non-planar curve without a
  /// given normal, and for a straight line).
  Vec3 plane_normal;
  /// Closed curves: +1 counter-clockwise, -1 clockwise, seen from the tip of plane_normal (by
  /// the sign of the enclosed area of the projection); 0 for open curves or without a plane.
  int orientation{0};
  double signed_area{0.0};  ///< enclosed area of the projection, > 0 counter-clockwise
};

/// Direction and sense of rotation. With a non-zero `plane_normal` the curve is projected along
/// it; otherwise the curve's own plane is used if it is planar within `tolerance`.
/// NurbsCurve::reversed() flips the direction.
[[nodiscard]] CurveDirection curve_direction(const NurbsCurve& c,
                                             const Vec3& plane_normal = {0, 0, 0},
                                             double tolerance = 1e-6);

// ------------------------------------------------------------------ continuity

struct ContinuityTolerance {
  double position{1e-6};   ///< largest gap of a joint (model units)
  double angle{1e-6};      ///< largest angle between the tangents for G1 (radians)
  double curvature{1e-4};  ///< largest difference of the curvature vectors for G2, relative to
                           ///< the larger curvature
};

/// Continuity of two curves at a joint (or of one curve at a parameter).
struct CurveContinuity {
  /// Geometric continuity: -1 the ends do not meet, 0 position (G0), 1 tangent (G1),
  /// 2 curvature (G2).
  int geometric{-1};
  /// Parametric continuity in the direction of travel through the joint: -1, 0 (C0),
  /// 1 (C1: equal first derivatives), 2 (C2: equal first and second derivatives).
  int parametric{-1};
  double gap{0.0};                            ///< distance of the two points
  double tangent_angle{0.0};                  ///< angle between the unit tangents (radians)
  double curvature_a{0.0}, curvature_b{0.0};  ///< curvatures on both sides
  double curvature_difference{0.0};           ///< |curvature vector a - curvature vector b|
};

/// Continuity at the joint of an end of a with an end of b. The curves are taken to run
/// through the joint from a to b: at End / Start as they are, otherwise reversed.
[[nodiscard]] CurveContinuity continuity(const NurbsCurve& a, CurveEnd end_a, const NurbsCurve& b,
                                         CurveEnd end_b, const ContinuityTolerance& tolerance = {});

/// Continuity at the pair of ends of a and b that are closest to each other.
[[nodiscard]] CurveContinuity continuity(const NurbsCurve& a, const NurbsCurve& b,
                                         const ContinuityTolerance& tolerance = {});

/// A parameter of a curve with its continuity there.
struct CurveBreak {
  double parameter{0.0};
  CurveContinuity continuity;
};

/// The interior knots of a curve where it is not curvature continuous (G2): corners, tangent
/// joints of lines and arcs, joints of joined curves, with the continuity found there from the
/// one-sided derivatives; for a closed curve also the seam, reported at domain_end. Empty for
/// a curve that is G2 throughout (a NURBS circle is, although it is only C1 at its knots).
[[nodiscard]] std::vector<CurveBreak> continuity_breaks(const NurbsCurve& c,
                                                        const ContinuityTolerance& tolerance = {});

// ------------------------------------------------------------------ duplicates

/// Groups of geometrically identical curves (SelDup): curves with the same shape within the
/// tolerance, whatever their direction, parametrisation, degree or knots. Every group lists
/// the indices of at least two curves in increasing order; the first is the one to keep.
/// Curves are compared by their lengths and by the distances of points of each to the other.
[[nodiscard]] std::vector<std::vector<std::size_t>> find_duplicate_curves(
    const std::vector<NurbsCurve>& curves, double tolerance = 1e-6);

}  // namespace nf
