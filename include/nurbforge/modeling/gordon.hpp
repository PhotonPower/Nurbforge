#pragma once
/// @file gordon.hpp
/// Coons patches and Gordon surfaces (Boolean sums, Piegl & Tiller section 10.5).

#include <vector>

#include "nurbforge/nurbs/curve.hpp"
#include "nurbforge/nurbs/surface.hpp"

namespace nf {

/// Bilinearly blended Coons patch from four boundary curves:
/// S(u, 0) = c0(u), S(u, 1) = c1(u), S(0, v) = d0(v), S(1, v) = d1(v) on [0, 1]^2.
/// Curve orientations are corrected automatically (d0 must connect the starts of c0 and c1,
/// d1 their ends); corners must agree within `tolerance`. Non-rational boundaries are
/// interpolated exactly. Rational boundaries (arcs, conics) are exact too when their weights
/// can be made to agree at the corners by scaling each curve's weights (e.g. arcs from
/// make_arc, whose end weights are 1): the Boolean sum is then taken in homogeneous
/// coordinates and the result is rational. Otherwise they are replaced by non-rational
/// approximations within `tolerance` first.
[[nodiscard]] NurbsSurface coons(const NurbsCurve& c0, const NurbsCurve& c1, const NurbsCurve& d0,
                                 const NurbsCurve& d1, double tolerance = 1e-7);

/// Gordon surface through a curve network: u_curves run in u and are ordered in v, v_curves run
/// in v and are ordered in u (at least two of each). Every u-curve must meet every v-curve,
/// with consistent parameters: the intersection parameters u_l (on all u-curves) and v_k (on
/// all v-curves) are determined from the network and checked within `tolerance`. The result
/// interpolates all network curves (domain [0, 1]^2); degrees across the curves are
/// min(3, count - 1). Rational curves are approximated as in coons().
[[nodiscard]] NurbsSurface gordon(const std::vector<NurbsCurve>& u_curves,
                                  const std::vector<NurbsCurve>& v_curves, double tolerance = 1e-7);

}  // namespace nf
