#pragma once
/// @file surface_tools.hpp
/// Convenience constructions on top of the core toolbox: rebuilding curves and surfaces with a
/// prescribed control net, tapered extrusions, ribbons, extrusions to a point and three-sided
/// patches.

#include <cstddef>
#include <vector>

#include "nurbforge/core/vec.hpp"
#include "nurbforge/nurbs/curve.hpp"
#include "nurbforge/nurbs/surface.hpp"

namespace nf {

/// Result of rebuild_curve.
struct RebuiltCurve {
  NurbsCurve curve;
  double max_error{0.0};  ///< max |c(t) - rebuilt(t)| at dense samples (same parameter)
};

/// Result of rebuild_surface.
struct RebuiltSurface {
  NurbsSurface surface;
  double max_error{0.0};  ///< max |S(u, v) - rebuilt(u, v)| on a dense grid (same parameters)
};

/// Rebuilds a curve as a non-rational B-spline of the given degree with `num_ctrl` control
/// points and uniform knots on the same domain: least-squares fit of dense samples at their
/// parameters (end points interpolated), so the parametrisation is kept. The reported error
/// compares points at equal parameters, an upper bound of the geometric deviation. Throws
/// GeometryError if num_ctrl < degree + 1.
[[nodiscard]] RebuiltCurve rebuild_curve(const NurbsCurve& curve, std::size_t num_ctrl,
                                         std::size_t degree = 3);

/// Rebuilds a surface as a non-rational B-spline surface with num_u x num_v control points,
/// degrees degree_u x degree_v and uniform knots on the same domain (least squares on a dense
/// grid, rows then columns; corners interpolated). The error is measured as for rebuild_curve on
/// a grid four times denser than the fit grid.
[[nodiscard]] RebuiltSurface rebuild_surface(const NurbsSurface& surface, std::size_t num_u,
                                             std::size_t num_v, std::size_t degree_u = 3,
                                             std::size_t degree_v = 3);

/// Tapered (drafted) extrusion of a planar profile along `direction` (normal to the profile's
/// plane; its length is the height): the top curve is the profile offset in its plane by
/// height * tan(angle), outwards for angle > 0 (the section grows), inwards for angle < 0, and
/// moved by `direction`. For an open profile "outwards" is the left side seen from the tip of
/// `direction`. The profile is split at its corners (tangent discontinuities); the offset
/// pieces are mitred (extended or cut where neighbours meet), and every piece gives one ruled
/// surface from the profile piece to its top piece. Offsets are fitted to `tolerance`
/// (offset_curve); line pieces give exact planar faces. |angle| must be below 90 degrees.
[[nodiscard]] std::vector<NurbsSurface> extrude_tapered(const NurbsCurve& profile,
                                                        const Vec3& direction, double angle,
                                                        double tolerance = 1e-7);

/// Planar ribbon along a curve: ruled surfaces between the curve and its offset by `width` in
/// the plane with normal `plane_normal` (width > 0: to the left, seen from the tip of the
/// normal). Corners are mitred as in extrude_tapered; one surface per smooth piece.
[[nodiscard]] std::vector<NurbsSurface> ribbon(const NurbsCurve& curve, double width,
                                               const Vec3& plane_normal, double tolerance = 1e-7);

/// Extrusion of a curve to a point (cone over the curve): the ruled surface from the curve
/// (u, as given; v = 0) to the apex (v = 1, a degenerate side with a limit normal). Exact, also
/// for rational curves (a circle gives an exact circular cone).
[[nodiscard]] NurbsSurface extrude_to_point(const NurbsCurve& curve, const Vec3& apex);

/// Three-sided patch bounded by c0, c1, c2 forming a closed loop (in any orientation; corners
/// must agree within `tolerance`): a Coons patch with c0 as one side and the corner opposite
/// to it as a degenerate side (v = 1), so the patch is regular except at that corner, where the
/// limit normal is used. Exact for non-rational curves; rational curves as in coons().
[[nodiscard]] NurbsSurface triangle_patch(const NurbsCurve& c0, const NurbsCurve& c1,
                                          const NurbsCurve& c2, double tolerance = 1e-7);

}  // namespace nf
