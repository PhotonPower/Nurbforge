#pragma once
/// @file trimmed_surface.hpp
/// Trimmed surfaces: a NURBS surface restricted by closed loops in its parameter domain.

#include <array>
#include <functional>
#include <vector>

#include "nurbforge/core/tolerance.hpp"
#include "nurbforge/nurbs/curve.hpp"
#include "nurbforge/nurbs/surface.hpp"

namespace nf {

/// Closed loop of curves in the parameter domain of a surface: x = u, y = v, z = 0. The curves
/// are consecutive (end of one = start of the next, last ends at the start of the first).
struct TrimLoop {
  std::vector<NurbsCurve> curves;

  /// Dense polyline approximation (closed: first point not repeated).
  [[nodiscard]] std::vector<std::array<double, 2>> polyline(std::size_t samples_per_span = 16) const;
  /// Signed area in the parameter plane (> 0 for counter-clockwise loops).
  [[nodiscard]] double signed_area() const;
  /// The same loop traversed in the opposite direction.
  [[nodiscard]] TrimLoop reversed() const;
};

/// Rectangle [u0, u1] x [v0, v1] as a counter-clockwise loop of four lines.
[[nodiscard]] TrimLoop rectangle_loop(double u0, double u1, double v0, double v1);
/// Exact circle of radius r about (cu, cv) in the parameter plane (counter-clockwise).
[[nodiscard]] TrimLoop circle_loop(double cu, double cv, double r);
/// Closed polygon through the given parameter points.
[[nodiscard]] TrimLoop polygon_loop(const std::vector<std::array<double, 2>>& points);

/// Surface restricted to the region inside `outer` and outside all `holes`. Loops are oriented
/// on construction: outer counter-clockwise, holes clockwise (in the (u, v) plane). Loops must be
/// closed (within tol.length in parameter units) and lie in the parameter domain.
class TrimmedSurface {
 public:
  /// Untrimmed: the outer loop is the parameter domain.
  explicit TrimmedSurface(NurbsSurface surface);
  TrimmedSurface(NurbsSurface surface, TrimLoop outer, std::vector<TrimLoop> holes = {},
                 const Tolerance& tol = {});

  [[nodiscard]] const NurbsSurface& surface() const noexcept { return surface_; }
  [[nodiscard]] const TrimLoop& outer() const noexcept { return outer_; }
  [[nodiscard]] const std::vector<TrimLoop>& holes() const noexcept { return holes_; }

  /// True if the region is the whole parameter domain (no holes, outer loop on the domain
  /// boundary).
  [[nodiscard]] bool is_untrimmed() const;

  /// True if (u, v) lies inside the trimmed region (on the boundary: either answer).
  [[nodiscard]] bool contains(double u, double v) const;

  /// Area of the trimmed region in 3D, by Green's theorem in the parameter plane:
  /// with F(u, v) = int |S_u x S_v| du, area = sum of the loop integrals of F dv (adaptive
  /// Gauss-Kronrod per knot span; holes, being clockwise, subtract).
  [[nodiscard]] double area(double rel_tol = 1e-9) const;

  /// Integral of f(u, v) over the trimmed region of the parameter plane (Green's theorem as for
  /// area(); include the area element |S_u x S_v| in f for surface integrals). `abs_tol` is the
  /// absolute error target of the loop integrals (needed when f nearly cancels).
  [[nodiscard]] double integrate(const std::function<double(double, double)>& f, double rel_tol,
                                 double abs_tol) const;

  /// The same surface with opposite orientation (u reversed: S'(u, v) = S(a + b - u, v)); the
  /// loops are mirrored accordingly. The normal S_u x S_v flips.
  [[nodiscard]] TrimmedSurface flipped() const;

  /// The boundary loops mapped to 3D: each trim curve composed with the surface, fitted to
  /// `tolerance` (cubic interpolation with adaptive refinement). Outer loop first.
  [[nodiscard]] std::vector<std::vector<NurbsCurve>> boundary_curves(double tolerance = 1e-7) const;

  /// Adds a hole (must lie inside the region).
  [[nodiscard]] TrimmedSurface with_hole(TrimLoop hole) const;

 private:
  NurbsSurface surface_;
  TrimLoop outer_;
  std::vector<TrimLoop> holes_;
  std::vector<std::vector<std::array<double, 2>>> poly_;  // outer + holes polylines (cached)
  void build_polylines();
};

/// Splits `surface` along its intersection curves with `cutter` into trimmed regions (the
/// p-curves on `surface` of surface/surface intersection). Open curves run from boundary to
/// boundary of the parameter domain (seams count as boundary) and must not cross each other;
/// closed curves become holes of the region around them and outer loops of new regions (nested
/// loops are supported). Returns the untrimmed surface if there is no intersection.
[[nodiscard]] std::vector<TrimmedSurface> split(const NurbsSurface& surface, const NurbsSurface& cutter,
                                                const Tolerance& tol = {});

/// Splits a trimmed face along parameter-space curves (general planar arrangement): the cuts
/// are clipped to the region, all curves are split at their mutual intersections, dangling
/// pieces are removed, and the faces of the resulting graph are traced (next half-edge clockwise
/// at every vertex). Counter-clockwise cycles inside the face become regions; clockwise cycles
/// (component boundaries) become holes of the smallest region containing them. Cuts may end on
/// the trim loops or on other cuts and may cross each other.
[[nodiscard]] std::vector<TrimmedSurface> split(const TrimmedSurface& face, const std::vector<NurbsCurve>& cuts,
                                                const Tolerance& tol = {});

/// Splits `surface` into regions along the given parameter-space curves (same rules as above).
[[nodiscard]] std::vector<TrimmedSurface> split(const NurbsSurface& surface,
                                                const std::vector<NurbsCurve>& pcurves,
                                                const Tolerance& tol = {});

}  // namespace nf
