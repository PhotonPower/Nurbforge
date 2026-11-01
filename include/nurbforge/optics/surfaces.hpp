#pragma once
/// @file surfaces.hpp
/// Optical surfaces (sag z as a function of the lateral position) as NURBS: conics (exact),
/// even/odd aspheres, Zernike and XY polynomial surfaces (fitted to a tolerance).
///
/// Conventions: the optical axis is z, the vertex lies at the origin, sag z(x, y) is measured
/// along +z. Curvature c = 1/R (R > 0: centre of curvature at +z), conic constant k
/// (k = 0 sphere, -1 paraboloid, < -1 hyperboloid, (-1, 0) prolate and > 0 oblate ellipsoid).

#include <cstddef>
#include <functional>
#include <vector>

#include "nurbforge/nurbs/curve.hpp"
#include "nurbforge/nurbs/surface.hpp"

namespace nf {

/// Rotationally symmetric sag
///   z(r) = c r^2 / (1 + sqrt(1 - (1+k) c^2 r^2)) + sum_i even[i] r^(2i+4) + sum_i odd[i] r^(i+1).
struct Asphere {
  double curvature{0.0};
  double conic{0.0};
  std::vector<double> even;  ///< A4, A6, A8, ... (coefficients of r^4, r^6, ...)
  std::vector<double> odd;   ///< B1, B2, B3, ... (coefficients of r^1, r^2, r^3, ...)

  [[nodiscard]] double sag(double r) const;
  /// Largest r for which the conic term is defined (infinite for k <= -1 or c = 0).
  [[nodiscard]] double max_radius() const;
  [[nodiscard]] bool is_conic() const { return even.empty() && odd.empty(); }
};

enum class ZernikeIndexing {
  Fringe,  ///< Fringe / University of Arizona ordering, unnormalised (Z4 = 2 rho^2 - 1)
  Noll,    ///< Noll ordering, normalised to unit RMS over the unit disk
};

/// Zernike polynomial Z_j(rho, theta) (j >= 1).
[[nodiscard]] double zernike(std::size_t j, ZernikeIndexing indexing, double rho, double theta);

/// Radial order n and azimuthal frequency m (m < 0: sine term) of Z_j.
void zernike_nm(std::size_t j, ZernikeIndexing indexing, int& n, int& m);

/// Zernike sag: base asphere plus sum_j coeffs[j-1] Z_j(r / norm_radius, theta).
struct ZernikeSag {
  Asphere base;
  double norm_radius{1.0};
  ZernikeIndexing indexing{ZernikeIndexing::Fringe};
  std::vector<double> coeffs;  ///< coefficient of Z_1, Z_2, ...
  [[nodiscard]] double sag(double x, double y) const;
};

/// XY polynomial sag: base asphere plus sum c * (x/norm_radius)^i (y/norm_radius)^j.
struct XYPolynomialSag {
  struct Term {
    int i{0}, j{0};
    double c{0.0};
  };
  Asphere base;
  double norm_radius{1.0};
  std::vector<Term> terms;
  [[nodiscard]] double sag(double x, double y) const;
};

/// Profile curve (r, 0, z(r)) for 0 <= r <= r_max in the xz-plane. Pure conics are exact
/// (piecewise rational quadratic); otherwise a cubic spline graph whose sag error is below
/// `tolerance` at dense check points (x(t) = r exactly, slope 0 at the vertex).
[[nodiscard]] NurbsCurve asphere_profile(const Asphere& a, double r_max, double tolerance = 1e-7);

/// Surface of revolution of asphere_profile about the z axis (u: azimuth, v: radius).
/// Exact for pure conics (sphere, ellipsoid, paraboloid, hyperboloid).
[[nodiscard]] NurbsSurface asphere_surface(const Asphere& a, double r_max, double tolerance = 1e-7);

/// Graph surface (x, y, f(x, y)) over [x0, x1] x [y0, y1]: tensor-product spline of the given
/// degree with x, y exactly linear; knots are inserted adaptively until |S_z - f| <= tolerance at
/// 8 x 8 check points per knot rectangle. Polynomials of degree <= `degree` are reproduced
/// exactly. Throws GeometryError if f is not finite somewhere or max_spans is exceeded.
[[nodiscard]] NurbsSurface sag_surface(const std::function<double(double, double)>& f, double x0,
                                       double x1, double y0, double y1, std::size_t degree = 3,
                                       double tolerance = 1e-7, std::size_t max_spans = 256);

/// Zernike surface over the square [-half_width, half_width]^2 (default: norm_radius).
[[nodiscard]] NurbsSurface zernike_surface(const ZernikeSag& z, double half_width = 0.0,
                                           double tolerance = 1e-7);

/// XY polynomial surface over the square [-half_width, half_width]^2 (default: norm_radius).
[[nodiscard]] NurbsSurface xy_polynomial_surface(const XYPolynomialSag& p, double half_width = 0.0,
                                                 double tolerance = 1e-7);

}  // namespace nf
