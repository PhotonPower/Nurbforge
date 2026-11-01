#include "nurbforge/modeling/loft.hpp"

#include <algorithm>

#include "nurbforge/core/error.hpp"
#include "nurbforge/modeling/construct.hpp"

namespace nf {

namespace {

std::vector<NurbsCurve> compatible(const std::vector<NurbsCurve>& sections) {
  if (sections.size() < 2) throw GeometryError("loft: need at least two sections");
  std::vector<NurbsCurve> c = sections;
  make_compatible(c);
  return c;
}

// Column i of the compatible sections: control points P_i of every section.
std::vector<std::vector<Vec3>> control_columns(const std::vector<NurbsCurve>& c) {
  const std::size_t n = c.front().num_control_points();
  std::vector<std::vector<Vec3>> cols(n, std::vector<Vec3>(c.size()));
  for (std::size_t k = 0; k < c.size(); ++k)
    for (std::size_t i = 0; i < n; ++i) cols[i][k] = c[k].control_points()[i];
  return cols;
}

NurbsSurface skin(const std::vector<NurbsCurve>& c, std::size_t degree_v,
                  const std::vector<double>& v) {
  const std::size_t K = c.size(), n = c.front().num_control_points();
  if (v.size() != K) throw GeometryError("loft: one parameter per section required");
  const std::size_t q = std::max<std::size_t>(1, std::min(degree_v, K - 1));
  // interpolate every column in homogeneous space: w*P and w with the same parameters/knots
  std::vector<Vec3> P;
  std::vector<double> W;
  KnotVector V;
  std::size_t nv = 0;
  std::vector<std::vector<Vec4>> rows(n);
  for (std::size_t i = 0; i < n; ++i) {
    std::vector<Vec3> wp(K), w(K);
    for (std::size_t k = 0; k < K; ++k) {
      const double wk = c[k].weights()[i];
      wp[k] = wk * c[k].control_points()[i];
      w[k] = {wk, 0.0, 0.0};
    }
    const NurbsCurve a = interpolate_curve(wp, q, v);
    const NurbsCurve b = interpolate_curve(w, q, v);
    if (i == 0) {
      V = a.knots();
      nv = a.num_control_points();
    }
    for (std::size_t j = 0; j < nv; ++j) {
      const double wj = b.control_points()[j].x;
      if (!(wj > 0.0))
        throw GeometryError("loft: interpolated weights are not positive (rational sections with "
                            "very different weights); reparametrise or use non-rational sections");
      rows[i].push_back(Vec4::from_point(a.control_points()[j] / wj, wj));
    }
  }
  for (std::size_t i = 0; i < n; ++i)
    for (std::size_t j = 0; j < nv; ++j) {
      P.push_back(rows[i][j].project());
      W.push_back(rows[i][j].w);
    }
  const NurbsCurve& c0 = c.front();
  return {c0.degree(), q, c0.knots(), std::move(V), n, nv, std::move(P), std::move(W)};
}

}  // namespace

std::vector<double> loft_parameters(const std::vector<NurbsCurve>& sections, Parametrization method) {
  return data_parameters(control_columns(compatible(sections)), method);
}

NurbsSurface loft(const std::vector<NurbsCurve>& sections, std::size_t degree_v,
                  Parametrization method) {
  const auto c = compatible(sections);
  return skin(c, degree_v, data_parameters(control_columns(c), method));
}

NurbsSurface loft(const std::vector<NurbsCurve>& sections, std::size_t degree_v,
                  const std::vector<double>& v_params) {
  return skin(compatible(sections), degree_v, v_params);
}

}  // namespace nf
