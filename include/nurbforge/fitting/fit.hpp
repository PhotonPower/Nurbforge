#pragma once
/// @file fit.hpp
/// Interpolation and least-squares approximation of point data (Piegl & Tiller, chapter 9).

#include <cstddef>
#include <optional>
#include <utility>
#include <vector>

#include "nurbforge/core/vec.hpp"
#include "nurbforge/geometry/plane.hpp"
#include "nurbforge/nurbs/curve.hpp"
#include "nurbforge/nurbs/surface.hpp"

namespace nf {

/// Choice of data parameters on [0, 1] (P&T Eq. 9.3–9.6).
enum class Parametrization {
  Uniform,      ///< equally spaced
  ChordLength,  ///< proportional to the distance between consecutive points
  Centripetal,  ///< proportional to the square root of the distance (robust at sharp turns)
};

/// Data parameters for the points (first 0, last 1). Throws GeometryError for fewer than two
/// points or coincident consecutive points (non-uniform parametrisations).
[[nodiscard]] std::vector<double> data_parameters(const std::vector<Vec3>& points,
                                                  Parametrization method);

/// Parameters on [0, 1] averaged over several point sequences of equal length (P&T Eq. 9.18,
/// 10.8); degenerate sequences (all points equal) are skipped.
[[nodiscard]] std::vector<double> data_parameters(const std::vector<std::vector<Vec3>>& sequences,
                                                  Parametrization method);

/// Global interpolation with explicit, strictly increasing parameters (domain
/// [params.front(), params.back()], knots by averaging).
[[nodiscard]] NurbsCurve interpolate_curve(const std::vector<Vec3>& points, std::size_t degree,
                                           const std::vector<double>& params);

/// Global interpolation (A9.1): non-rational B-spline of the given degree on [0, 1] through all
/// points, knots by averaging (Eq. 9.8). Optional end derivatives {C'(0), C'(1)} (w.r.t. the
/// parameter on [0, 1]) add two control points (Eq. 9.21/9.22). Needs at least degree+1 points.
[[nodiscard]] NurbsCurve interpolate_curve(
    const std::vector<Vec3>& points, std::size_t degree,
    Parametrization method = Parametrization::Centripetal,
    const std::optional<std::pair<Vec3, Vec3>>& end_derivatives = std::nullopt);

/// Least-squares approximation (Eq. 9.63–9.69) with n_ctrl control points; the first and last
/// point are interpolated. degree+1 <= n_ctrl <= number of points.
[[nodiscard]] NurbsCurve approximate_curve(const std::vector<Vec3>& points, std::size_t degree,
                                           std::size_t n_ctrl,
                                           Parametrization method = Parametrization::Centripetal);

/// Least-squares approximation with given knot vector (clamped, degree+1 <= n_ctrl <= number of
/// points) and data parameters (inside the knot domain, first/last at the domain ends); the
/// first and last point are interpolated. Reproduces data from the spline space exactly, e.g.
/// x(t) = t stays exactly linear (used for graph fits).
[[nodiscard]] NurbsCurve approximate_curve(const std::vector<Vec3>& points, std::size_t degree,
                                           const KnotVector& knots,
                                           const std::vector<double>& params);

/// Least-squares approximation with the smallest number of control points such that every data
/// point lies within `tolerance` of the curve point at its parameter (an upper bound of the
/// point-to-curve distance). Falls back to interpolation if necessary.
[[nodiscard]] NurbsCurve approximate_curve_tol(
    const std::vector<Vec3>& points, std::size_t degree, double tolerance,
    Parametrization method = Parametrization::Centripetal);

/// Global surface interpolation (A9.4) through a grid of points, row-major
/// points[i * num_v + j] (i along u). Parameters are averaged over rows/columns; degenerate rows
/// (e.g. a pole, all points equal) are skipped in the averaging.
[[nodiscard]] NurbsSurface interpolate_surface(
    const std::vector<Vec3>& points, std::size_t num_u, std::size_t num_v, std::size_t degree_u,
    std::size_t degree_v, Parametrization method = Parametrization::Centripetal);

/// Least-squares surface approximation of a point grid with ctrl_u x ctrl_v control points
/// (curve approximation through all rows, then through the resulting columns; corners are
/// interpolated).
[[nodiscard]] NurbsSurface approximate_surface(
    const std::vector<Vec3>& points, std::size_t num_u, std::size_t num_v, std::size_t degree_u,
    std::size_t degree_v, std::size_t ctrl_u, std::size_t ctrl_v,
    Parametrization method = Parametrization::Centripetal);

/// Best-fit line through points (orthogonal regression: minimises the sum of squared
/// perpendicular distances).
struct LineFit {
  Vec3 point;                 ///< centroid of the points (on the line)
  Vec3 direction;             ///< unit direction; sign: its largest-magnitude component is positive
  double max_deviation{0.0};  ///< largest perpendicular distance of a point from the line
  double rms{0.0};            ///< root mean square of the perpendicular distances
};

/// Best-fit plane through points (orthogonal regression).
struct PlaneFit {
  Plane plane;  ///< through the centroid; unit normal, largest-magnitude component positive
  double max_deviation{0.0};  ///< largest distance of a point from the plane
  double rms{0.0};            ///< root mean square of the distances
};

/// Line through points by orthogonal regression: the centroid and the eigenvector of the largest
/// eigenvalue of the scatter matrix (symmetric Jacobi rotations). Needs at least 2 distinct
/// points; throws GeometryError otherwise.
[[nodiscard]] LineFit fit_line(const std::vector<Vec3>& points);

/// Plane through points by orthogonal regression: the centroid and the eigenvector of the
/// smallest eigenvalue of the scatter matrix. Needs at least 3 points that are not collinear;
/// throws GeometryError otherwise.
[[nodiscard]] PlaneFit fit_plane(const std::vector<Vec3>& points);

}  // namespace nf
