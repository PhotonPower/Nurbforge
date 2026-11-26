#pragma once
/// @file unroll.hpp
/// Developable surfaces: a developability check and the unrolling into the plane (flat
/// patterns for sheet metal, foils and flexible substrates).

#include "nurbforge/nurbs/surface.hpp"

namespace nf {

/// Parameter direction of the straight lines of a ruled surface.
enum class RulingDirection {
  None,  ///< not ruled along a parameter direction
  U,     ///< the lines v = const are straight (the surface is linear in u)
  V,     ///< the lines u = const are straight (the surface is linear in v)
};

/// Result of developability().
struct Developability {
  /// Ruled along a parameter direction: its parameter lines there are straight and linearly
  /// parametrised (extrusions, ruled surfaces, cones to a point, cylinders and cones of
  /// revolution).
  RulingDirection rulings{RulingDirection::None};
  /// Ruled, and the normal is the same all along every ruling (within the angular tolerance):
  /// the surface can be unrolled into the plane without stretching.
  bool developable{false};
  /// Largest angle (radians) between the normals at two points of one ruling, at the sampled
  /// rulings; 0 for a developable surface. For a surface that is not ruled: 0.
  double max_twist{0.0};
  /// Largest |Gaussian curvature| at the sample points (0 for a developable surface).
  double max_gaussian_curvature{0.0};
};

/// Checks whether a surface is developable. Rulings are searched along the parameter
/// directions only (v first): the surface must agree with the straight lines between its two
/// boundary curves there within 1e-9 of its size. Normals are compared at 9 points of each of
/// 16 rulings per knot span.
[[nodiscard]] Developability developability(const NurbsSurface& surface,
                                            double angular_tolerance = 1e-6);

/// Result of unroll().
struct UnrolledSurface {
  /// The flat pattern in the plane z = 0 on the parameter domain of the surface:
  /// pattern(u, v) is the image of S(u, v). Its normal S_u x S_v is +z, so it is seen from the
  /// side of the surface normal (no mirror image). The directrix (the middle of the rulings)
  /// starts at the origin in the direction +x.
  NurbsSurface pattern;
  /// Largest deviation of the pattern from the exact development at the check points of the
  /// fit (model units, <= tolerance).
  double max_error{0.0};
};

/// Unrolls a developable surface into the plane without stretching: lengths, angles and areas
/// on the surface are kept, so curves on it map to the pattern through their parameters,
/// pattern(u(t), v(t)).
///
/// The curve through the middle of the rulings is developed by its arc length and its geodesic
/// curvature (both integrated by Gauss-Legendre quadrature to the tolerance), and the rulings
/// are laid out at the angles they make with it. The pattern is linear along the rulings and
/// a cubic spline across them, interpolated until it is within `tolerance` of the exact
/// development, piecewise between the knots of the surface; it keeps the parametrisation of
/// the surface, so a cylinder of revolution gives a rectangle and a cone a circular sector
/// within the tolerance (the parametrisation of a rational circle by its arc length is not
/// rational).
///
/// Throws GeometryError if the surface is not ruled along a parameter direction, if it is not
/// developable within `angular_tolerance` (see developability()), or if the tolerance cannot
/// be reached. Overlaps of the pattern (a surface winding more than once) are not detected.
[[nodiscard]] UnrolledSurface unroll(const NurbsSurface& surface, double tolerance = 1e-7,
                                     double angular_tolerance = 1e-6);

}  // namespace nf
