#pragma once
/// @file drawing.hpp
/// 2D drawings of 3D models: visibility of curves against shells (hidden-line removal), the
/// visible and hidden lines of shells projected onto a plane (Make2D), and the outline of shells
/// seen along a direction.

#include <cstddef>
#include <vector>

#include "nurbforge/analysis/silhouette.hpp"
#include "nurbforge/core/vec.hpp"
#include "nurbforge/nurbs/curve.hpp"
#include "nurbforge/topology/region.hpp"
#include "nurbforge/topology/shell.hpp"

namespace nf {

struct VisibilityOptions {
  /// Model-space tolerance: a surface hides a point only if the point lies more than this off
  /// the surface's tangent plane where the view ray meets it. Curves lying on an occluding
  /// surface (the edges of a face, silhouettes) are therefore not hidden by that surface.
  double tolerance{1e-6};
  /// Visibility is sampled at max(samples, 4 x knot spans) + 1 points per curve; the parameters
  /// where it changes are then refined by bisection. A hidden or visible stretch shorter than
  /// the sample spacing can be missed.
  std::size_t samples{64};
};

/// A piece of a curve with constant visibility.
struct CurvePiece {
  NurbsCurve curve;  ///< the part of the source curve on [t0, t1] (same parameter)
  bool visible{true};
  std::size_t source{0};    ///< index of the source curve
  double t0{0.0}, t1{0.0};  ///< parameter interval on the source curve
};

/// Visibility of curves against shells: every curve is split at the parameters where its
/// visibility changes, and each piece is classified as visible or hidden.
///
/// A point is hidden if the ray from it towards the viewer (parallel view: against the view
/// direction; perspective: to the eye point) crosses a face of an occluder, see
/// VisibilityOptions::tolerance for what counts as crossing; rays grazing a face or lying in
/// it do not hide. Pieces are returned in the order of the curves and of the parameter.
[[nodiscard]] std::vector<CurvePiece> visibility(const std::vector<NurbsCurve>& curves,
                                                 const View& view,
                                                 const std::vector<Shell>& occluders,
                                                 const VisibilityOptions& options = {});

struct DrawingOptions {
  double tolerance{1e-6};     ///< fit tolerance of silhouettes, visibility tolerance (model units)
  std::size_t samples{64};    ///< VisibilityOptions::samples
  std::size_t resolution{8};  ///< SilhouetteOptions::resolution
  /// First axis s of the drawing plane (projected into the plane and normalised); zero: chosen
  /// automatically. The second axis is t = normal x s, the normal points to the viewer.
  Vec3 s_axis{0, 0, 0};
  Vec3 origin{0, 0, 0};  ///< point mapped to (0, 0)
  /// Two projected curves closer than this are the same line of the drawing (duplicates are
  /// removed); also the tolerance for chaining the outline. 0: 10 x tolerance.
  double join_tolerance{0.0};
  bool hidden{true};  ///< make2d: also return the hidden lines
  /// Also draw smooth edges, where the two faces share their tangent plane (the seam of a
  /// cylinder, the borders of a fillet). They are not visible features of the shape.
  bool tangent_edges{false};
};

enum class DrawingCurveKind {
  Edge,        ///< edge of a shell
  Silhouette,  ///< silhouette curve of a face
};

/// One line of a drawing.
struct DrawingCurve {
  NurbsCurve curve;    ///< in the drawing plane: (s, t, 0)
  NurbsCurve curve3d;  ///< the same piece in 3D (same parameter)
  bool visible{true};
  DrawingCurveKind kind{DrawingCurveKind::Edge};
  std::size_t shell{0};  ///< index of the source shell
  std::size_t index{0};  ///< edge index (Shell::edges) or face index of the silhouette
};

/// Drawing of shells seen along a direction, in the plane coordinates
/// p = origin + s * s_axis + t * t_axis (seen from the viewer, s points right and t up).
struct Drawing {
  Vec3 origin, s_axis, t_axis, normal;
  std::vector<DrawingCurve> visible;
  std::vector<DrawingCurve> hidden;
};

/// Make2D: all edges and silhouettes of the shells for a parallel view along `direction` (from
/// the viewer into the scene), split by visibility and projected onto a plane normal to the
/// direction (an exact affine map of the control points).
///
/// Curves are first split where their projection folds back (a circle seen edge-on), pieces
/// projecting to a point are dropped, and duplicates are removed: a piece lying on a line
/// already drawn is cut out, visible lines taking precedence over hidden ones (the far edge of
/// a box under the near one, a silhouette on an edge).
[[nodiscard]] Drawing make2d(const std::vector<Shell>& shells, const Vec3& direction,
                             const DrawingOptions& options = {});

/// Outline of shells seen along a direction.
struct Outline {
  Vec3 origin, s_axis, t_axis, normal;  ///< plane coordinates as for Drawing
  /// Boundary of the union of the projections: closed loops of curves (s, t, 0) with the
  /// shells' projection on the left, i.e. outer contours counter-clockwise and holes clockwise.
  Region loops;
};

/// The contour of shells seen along `direction`: the boundary of the union of their parallel
/// projections onto a plane normal to the direction. The projected edges and silhouettes are
/// split at their mutual intersections; a piece belongs to the outline if the shells cover
/// exactly one of its two sides (view rays next to it on either side). The pieces are chained
/// into closed loops. Throws GeometryError if a chain does not close.
[[nodiscard]] Outline outline(const std::vector<Shell>& shells, const Vec3& direction,
                              const DrawingOptions& options = {});

}  // namespace nf
