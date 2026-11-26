#pragma once
// Internal: adaptive cubic interpolation of several curves sharing one parameter (not
// installed).

#include <algorithm>
#include <functional>
#include <map>
#include <string>
#include <vector>

#include "nurbforge/core/error.hpp"
#include "nurbforge/fitting/fit.hpp"
#include "nurbforge/nurbs/curve.hpp"

namespace nf::detail {

/// Cubic interpolation of a family of curves with a common parameter: `exact(t)` gives one
/// point per curve, `error(curves, t, exact(t))` the deviation of the interpolants at t. The
/// curves interpolate the exact points at the parameters T (so they have equal knot vectors);
/// the midpoints of the intervals where the error exceeds `tol` are added until none is left.
/// Throws GeometryError naming `who` if the tolerance is not reached.
inline std::vector<NurbsCurve> fit_family(
    std::vector<double> T, const std::function<std::vector<Vec3>(double)>& exact,
    const std::function<double(const std::vector<NurbsCurve>&, double, const std::vector<Vec3>&)>&
        error,
    double tol, const std::string& who) {
  std::map<double, std::vector<Vec3>> cache;
  const auto at = [&](double t) -> const std::vector<Vec3>& {
    auto it = cache.find(t);
    if (it == cache.end()) it = cache.emplace(t, exact(t)).first;
    return it->second;
  };
  const std::size_t ncurves = at(T.front()).size();
  for (int round = 0; round < 16; ++round) {
    std::vector<std::vector<Vec3>> pts(ncurves);
    for (const double t : T)
      for (std::size_t k = 0; k < ncurves; ++k) pts[k].push_back(at(t)[k]);
    const std::size_t deg = std::min<std::size_t>(3, T.size() - 1);
    std::vector<NurbsCurve> curves;
    for (std::size_t k = 0; k < ncurves; ++k) curves.push_back(interpolate_curve(pts[k], deg, T));
    std::vector<double> add;
    for (std::size_t i = 0; i + 1 < T.size(); ++i) {
      const double tm = 0.5 * (T[i] + T[i + 1]);
      if (error(curves, tm, at(tm)) > tol) add.push_back(tm);
    }
    if (add.empty()) return curves;
    if (T.size() + add.size() > 8193) break;
    T.insert(T.end(), add.begin(), add.end());
    std::sort(T.begin(), T.end());
  }
  throw GeometryError(who + ": tolerance " + std::to_string(tol) + " not reached");
}

}  // namespace nf::detail
