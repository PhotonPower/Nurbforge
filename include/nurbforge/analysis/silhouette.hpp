#pragma once
/// @file silhouette.hpp
/// Views (parallel or perspective) and silhouette curves of surfaces, trimmed faces and shells:
/// the curves where the surface turns away from the viewer.

#include <cstddef>
#include <optional>
#include <vector>

#include "nurbforge/core/vec.hpp"
#include "nurbforge/nurbs/curve.hpp"
#include "nurbforge/nurbs/surface.hpp"
#include "nurbforge/topology/shell.hpp"
#include "nurbforge/topology/trimmed_surface.hpp"

namespace nf {

/// A view of the model: parallel along `direction`, or perspective from the point `eye`.
struct View {
  /// Viewing direction of a parallel view, from the viewer into the scene (need not be
  /// normalised). Ignored if `eye` is set.
  Vec3 direction{0, 0, -1};
  /// Eye point of a perspective view.
  std::optional<Vec3> eye{};

  [[nodiscard]] static View parallel(const Vec3& direction);
  [[nodiscard]] static View perspective(const Vec3& eye);
  /// Unit vector from p towards the viewer. Throws GeometryError for a zero direction or if p
  /// is the eye point.
  [[nodiscard]] Vec3 to_viewer(const Vec3& p) const;
};

struct SilhouetteOptions {
  /// Accuracy of the curves (model units): the 3D curve and the image of the p-curve stay
  /// within this distance of exact silhouette points, verified between all fit points.
  double tolerance{1e-6};
  /// Cells per knot span of the surface (each direction) of the grid on which sign changes of
  /// N . w are searched. A closed silhouette smaller than a cell, or two folds closer than a
  /// cell, can be missed; raise it for surfaces with fine detail inside single knot spans.
  std::size_t resolution{8};
};

/// One silhouette curve. `curve` and `pcurve` share their parameter.
struct SilhouetteCurve {
  NurbsCurve curve;     ///< 3D curve on the surface
  NurbsCurve pcurve;    ///< (u, v, 0) in the parameter domain of the surface
  bool closed{false};   ///< the 3D curve is closed
  bool crease{false};   ///< a crease of the surface (C0 knot line) between a part facing the
                        ///< viewer and a part facing away, not a smooth fold
  std::size_t face{0};  ///< index of the source face (shells)
  /// Parameters where the view ray is tangent to the curve: its projection has a cusp there
  /// and the silhouette changes between the near and the far side of a fold.
  std::vector<double> cusps;
};

/// Silhouette curves of a surface for a view: the points where the surface normal N is
/// perpendicular to the view ray, N . d = 0 (parallel view) or N . (S - eye) = 0 (perspective).
///
/// The function g(u, v) = N . w (w the unit vector to the viewer, N the unit normal with its
/// limit at poles) is sampled on a grid over the parameter domain (options.resolution cells
/// per knot span); its zero set is traced through the cells with a sign change, the points
/// are moved onto the exact zero set by Newton steps and interpolated (cubic) until the curve
/// and the image of its p-curve are within options.tolerance of exact silhouette points.
/// Straight silhouettes (cylinders and cones seen from the side) are therefore straight to
/// round-off. Curves are split where they cross the boundary of the parameter domain (seams).
/// Knot lines across which the normal jumps (creases) are reported with `crease` set where the
/// two sides face opposite ways.
///
/// A surface seen edge-on everywhere (a planar face containing the view direction, a cylinder
/// seen along its axis) has no silhouette curves: its outline is its boundary. Isolated points
/// with N . w = 0 are not reported.
[[nodiscard]] std::vector<SilhouetteCurve> silhouette(const NurbsSurface& surface, const View& view,
                                                      const SilhouetteOptions& options = {});

/// Silhouette curves of a trimmed face: those of its surface, clipped to the trimmed region.
[[nodiscard]] std::vector<SilhouetteCurve> silhouette(const TrimmedSurface& face, const View& view,
                                                      const SilhouetteOptions& options = {});

/// Silhouette curves of all faces of a shell (SilhouetteCurve::face gives the face).
[[nodiscard]] std::vector<SilhouetteCurve> silhouette(const Shell& shell, const View& view,
                                                      const SilhouetteOptions& options = {});

}  // namespace nf
