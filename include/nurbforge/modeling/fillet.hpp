#pragma once
/// @file fillet.hpp
/// Constant-radius fillets: circular arcs between two planar curves and rolling-ball fillet
/// surfaces between two surfaces.
///
/// Both use the offset construction: the centre of a fillet of radius r lies at distance r from
/// both inputs, i.e. on the intersection of their offsets by r. The offsets (modeling/offset.hpp)
/// only provide the start values; centres and contact points are then solved exactly
/// (Newton on the foot-point conditions), so the fillet does not inherit the offset fitting error.

#include <cstddef>
#include <optional>

#include "nurbforge/core/vec.hpp"
#include "nurbforge/nurbs/curve.hpp"
#include "nurbforge/nurbs/surface.hpp"

namespace nf {

/// Options of fillet_curves.
struct CurveFilletOptions {
  /// Side of the fillet centre relative to curve a: +1 = left of a's direction of travel seen
  /// from the tip of plane_normal (the side of offset_curve(a, +r, plane_normal)), -1 = right,
  /// 0 = automatic (the corner formed by the kept parts of a and b, see fillet_curves).
  int side_a{0};
  int side_b{0};  ///< same for curve b (relative to b's direction of travel)
  /// Trim a and b at the tangent points (CurveFillet::trimmed_a / trimmed_b).
  bool trim{true};
  /// If the curves intersect more than once: use the intersection nearest to this point
  /// (default: nearest to the midpoint of a's end point and b's start point).
  std::optional<Vec3> corner{};
  /// Fit tolerance of the auxiliary offset curves (start values only; the result is exact).
  double tolerance{1e-7};
};

/// Result of fillet_curves.
struct CurveFillet {
  /// Exact rational arc of radius r (make_arc) from point_a to point_b, oriented so that
  /// trimmed_a, arc, trimmed_b form a G1 chain (join(join(trimmed_a, arc), trimmed_b)).
  NurbsCurve arc;
  Vec3 center;         ///< centre of the arc
  double param_a{0.0}; ///< parameter of the tangent point on a
  double param_b{0.0}; ///< parameter of the tangent point on b
  Vec3 point_a;        ///< a(param_a) = arc start
  Vec3 point_b;        ///< b(param_b) = arc end
  /// The kept part of a (the side of param_a away from the corner), oriented to end at point_a;
  /// reversed with respect to a if the kept part lies after param_a. Empty if options.trim is
  /// false.
  std::optional<NurbsCurve> trimmed_a;
  /// The kept part of b, oriented to start at point_b (reversed if it lies before param_b).
  std::optional<NurbsCurve> trimmed_b;
};

/// Circular fillet of the given radius between two coplanar, intersecting curves.
///
/// The corner X is an intersection of a and b (the one nearest to options.corner). On each
/// curve the **kept part** is the side of X with the longer parameter interval (for a polyline
/// corner where a ends and b starts at X: a before X, b after X). With automatic sides the
/// fillet lies in the sector spanned by the kept parts: its centre is to the side of a on which
/// b's kept part lies, and vice versa. The centre is the intersection of
/// offset_curve(a, side_a * r) and offset_curve(b, side_b * r) nearest to X, refined by Newton on
/// a(s) + d_a N_a(s) = b(t) + d_b N_b(t) (exact unit normals, N = normalize(n x c')), so the
/// tangent points are the exact foot points of the centre on both curves. Exact for lines
/// (tangent points at distance r / tan(theta / 2) from the corner for an angle theta between
/// the kept parts) and circles.
///
/// Throws GeometryError if the curves do not lie in a common plane normal to plane_normal, do
/// not intersect, are tangent at the corner, if r is not positive, if the offsets do not
/// intersect (fillet too large for the curves, or r beyond a radius of curvature on the concave
/// side), or if a tangent point falls on the far side of the corner (trimmed part empty).
[[nodiscard]] CurveFillet fillet_curves(const NurbsCurve& a, const NurbsCurve& b, double radius,
                                        const Vec3& plane_normal = {0.0, 0.0, 1.0},
                                        const CurveFilletOptions& options = {});

/// Options of fillet_surfaces.
struct SurfaceFilletOptions {
  /// Side of the ball centre relative to s1: centre = contact point + side1 * r * N1 with the
  /// unit normal N1 = normalize(S1_u x S1_v) (the side of offset_surface(s1, side1 * r)).
  /// Must be +1 or -1.
  int side1{1};
  int side2{1};  ///< same for s2
  /// Deviation of the fillet surface from the ideal rolling-ball surface (and of spine, contact
  /// curves and mapped p-curves from their exact counterparts), verified at check points
  /// between the sections.
  double tolerance{1e-6};
  /// Fit tolerance of the auxiliary offset surfaces (start values only); 0: tolerance / 10.
  double offset_tolerance{0.0};
  std::size_t max_sections{1025};
};

/// Result of fillet_surfaces. All curves share the parameter t of the spine.
struct SurfaceFillet {
  /// u in [0, 1]: circular arcs of radius r from the contact on s1 (u = 0) to the contact on s2
  /// (u = 1) around the spine point, in the plane normal to the spine tangent; v = t.
  NurbsSurface surface;
  NurbsCurve spine;     ///< path of the ball centre
  NurbsCurve contact1;  ///< contact curve on s1 (= surface boundary u = 0)
  NurbsCurve contact2;  ///< contact curve on s2 (= surface boundary u = 1)
  NurbsCurve pcurve1;   ///< (u1, v1, 0) of contact1 in the parameter domain of s1
  NurbsCurve pcurve2;   ///< (u2, v2, 0) of contact2 in the parameter domain of s2
};

/// Rolling-ball fillet of radius r between two surfaces.
///
/// Spine: the intersection curve of offset_surface(s1, side1 * r) and offset_surface(s2,
/// side2 * r) (the longest one if there are several), used as start values. At every section
/// parameter t the exact ball centre C(t) is solved in the plane through spine0(t) normal to the
/// spine tangent (Newton on |C - S1| = r, |C - S2| = r with closest-point projections); contact
/// points are the foot points (with parameters for the p-curves). Sections are the shorter
/// circular arcs (rational, exact) from contact 1 to contact 2 around C; they are skinned
/// (loft, cubic in v) and the section count is refined (interval bisection) until the deviation
/// from the ideal rolling-ball surface at 3 x 17 check points per interval is below
/// options.tolerance. Spine and p-curves are cubic interpolations of the exact section data,
/// contact curves are the boundary iso-curves of the fillet surface.
///
/// The fillet is untrimmed and extends over the whole spine (until the offsets leave a
/// parameter domain); s1 and s2 are not trimmed. Throws GeometryError if the offsets do not
/// intersect, the ball cannot be placed (Newton fails, e.g. at tangential configurations), or
/// the tolerance cannot be met with max_sections.
[[nodiscard]] SurfaceFillet fillet_surfaces(const NurbsSurface& s1, const NurbsSurface& s2,
                                            double radius,
                                            const SurfaceFilletOptions& options = {});

}  // namespace nf
