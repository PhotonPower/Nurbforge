#pragma once
/// @file loft.hpp
/// Lofting (skinning) a surface through section curves (Piegl & Tiller, section 10.3).

#include <cstddef>
#include <vector>

#include "nurbforge/fitting/fit.hpp"
#include "nurbforge/nurbs/curve.hpp"
#include "nurbforge/nurbs/surface.hpp"

namespace nf {

/// Parameters v_k of the sections (first 0, last 1), averaged over the columns of the
/// compatible control polygons (P&T Eq. 10.8). The sections are made compatible internally.
[[nodiscard]] std::vector<double> loft_parameters(const std::vector<NurbsCurve>& sections,
                                                  Parametrization method = Parametrization::Centripetal);

/// Skinned surface through the section curves, in the given order. u follows the sections
/// (domain [0, 1], made compatible: common degree and knots), v runs across them with degree
/// min(degree_v, number of sections - 1). S(u, v_k) = C_k(u) for the section parameters v_k of
/// loft_parameters. Rational sections are interpolated in homogeneous space; throws
/// GeometryError if that would produce non-positive weights, or for fewer than two sections.
[[nodiscard]] NurbsSurface loft(const std::vector<NurbsCurve>& sections, std::size_t degree_v = 3,
                                Parametrization method = Parametrization::Centripetal);

/// As above with explicit, strictly increasing section parameters (the v domain is
/// [v_params.front(), v_params.back()]).
[[nodiscard]] NurbsSurface loft(const std::vector<NurbsCurve>& sections, std::size_t degree_v,
                                const std::vector<double>& v_params);

}  // namespace nf
