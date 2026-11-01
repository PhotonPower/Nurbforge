#include "nurbforge/modeling/construct.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>

#include "nurbforge/core/error.hpp"
#include "nurbforge/geometry/primitives.hpp"

namespace nf {

NurbsSurface extrude(const NurbsCurve& profile, const Vec3& direction) {
  if (norm(direction) < tol::kLength) throw GeometryError("extrusion direction is zero");
  const std::size_t n = profile.num_control_points();
  std::vector<Vec3> P;
  std::vector<double> W;
  P.reserve(2 * n);
  W.reserve(2 * n);
  for (std::size_t i = 0; i < n; ++i) {
    P.push_back(profile.control_points()[i]);
    P.push_back(profile.control_points()[i] + direction);
    W.push_back(profile.weights()[i]);
    W.push_back(profile.weights()[i]);
  }
  return {profile.degree(), 1, profile.knots(), {0, 0, 1, 1}, n, 2, std::move(P), std::move(W)};
}

NurbsSurface revolve(const NurbsCurve& profile, const Vec3& S, const Vec3& axis_dir, double theta) {
  if (!(theta > 0.0) || theta > 2.0 * std::numbers::pi + tol::kAngular)
    throw GeometryError("revolution angle must be in (0, 2*pi]");
  const Vec3 T = normalized(axis_dir);
  const std::size_t narcs = detail::arc_segment_count(theta);
  const double dtheta = theta / static_cast<double>(narcs);
  const double wm = std::cos(dtheta / 2.0);
  const std::size_t nu = 2 * narcs + 1;
  const std::size_t nv = profile.num_control_points();

  std::vector<double> cosines(narcs + 1), sines(narcs + 1);
  for (std::size_t i = 1; i <= narcs; ++i) {
    cosines[i] = std::cos(dtheta * static_cast<double>(i));
    sines[i] = std::sin(dtheta * static_cast<double>(i));
  }

  std::vector<Vec3> P(nu * nv);
  std::vector<double> W(nu * nv);
  auto idx = [nv](std::size_t i, std::size_t j) { return i * nv + j; };

  for (std::size_t j = 0; j < nv; ++j) {
    const Vec3 Pj = profile.control_points()[j];
    const double wj = profile.weights()[j];
    const Vec3 O = S + dot(Pj - S, T) * T;  // projection onto axis
    Vec3 X = Pj - O;
    const double r = norm(X);
    if (r < tol::kLength) {  // point on the axis: degenerate (pole)
      for (std::size_t i = 0; i < nu; ++i) {
        P[idx(i, j)] = Pj;
        W[idx(i, j)] = (i % 2 == 1) ? wm * wj : wj;
      }
      continue;
    }
    X = X / r;
    const Vec3 Y = cross(T, X);
    P[idx(0, j)] = Pj;
    W[idx(0, j)] = wj;
    Vec3 P0 = Pj, T0 = Y;
    std::size_t index = 0;
    for (std::size_t i = 1; i <= narcs; ++i) {
      const Vec3 P2 = O + r * cosines[i] * X + r * sines[i] * Y;
      const Vec3 T2 = -sines[i] * X + cosines[i] * Y;
      P[idx(index + 2, j)] = P2;
      W[idx(index + 2, j)] = wj;
      P[idx(index + 1, j)] = detail::intersect_lines(P0, T0, P2, T2);
      W[idx(index + 1, j)] = wm * wj;
      index += 2;
      P0 = P2;
      T0 = T2;
    }
  }
  return {2, profile.degree(), detail::arc_knots(narcs), profile.knots(), nu, nv,
          std::move(P), std::move(W)};
}

NurbsCurve reparametrize(const NurbsCurve& c, double a, double b) {
  if (!(b > a)) throw GeometryError("reparametrize: need a < b");
  const double c0 = c.domain_start(), c1 = c.domain_end();
  KnotVector k = c.knots();
  for (double& x : k) x = a + (x - c0) * (b - a) / (c1 - c0);
  return {c.degree(), std::move(k), c.control_points(), c.weights()};
}

void make_compatible(NurbsCurve& a, NurbsCurve& b) {
  std::vector<NurbsCurve> v{a, b};
  make_compatible(v);
  a = std::move(v[0]);
  b = std::move(v[1]);
}

void make_compatible(std::vector<NurbsCurve>& curves) {
  if (curves.empty()) return;
  std::size_t p = 0;
  for (auto& c : curves) {
    c = reparametrize(c.clamped(), 0.0, 1.0);
    p = std::max(p, c.degree());
  }
  for (auto& c : curves)
    if (c.degree() < p) c = c.elevate_degree(p - c.degree());
  // union of all knot values with maximal multiplicity
  std::vector<double> vals;
  std::vector<std::size_t> mults;
  for (const auto& c : curves) {
    const auto [cv, cm] = unique_knots(c.knots());
    for (std::size_t i = 0; i < cv.size(); ++i) {
      std::size_t k = 0;
      while (k < vals.size() && std::abs(vals[k] - cv[i]) >= tol::kParametric) ++k;
      if (k == vals.size()) {
        vals.push_back(cv[i]);
        mults.push_back(0);
      }
      mults[k] = std::max(mults[k], static_cast<std::size_t>(cm[i]));
    }
  }
  for (auto& c : curves) {
    std::vector<double> X;
    for (std::size_t k = 0; k < vals.size(); ++k)
      for (auto have = knot_multiplicity(c.knots(), vals[k]); have < mults[k]; ++have) X.push_back(vals[k]);
    c = c.refine_knots(X);
  }
  // identical knot values (refinement snaps to existing knots, tiny differences may remain)
  for (auto& c : curves)
    c = NurbsCurve(c.degree(), curves.front().knots(), c.control_points(), c.weights());
}

NurbsSurface ruled(const NurbsCurve& c0, const NurbsCurve& c1) {
  NurbsCurve a = c0, b = c1;
  make_compatible(a, b);
  const std::size_t n = a.num_control_points();
  std::vector<Vec3> P;
  std::vector<double> W;
  for (std::size_t i = 0; i < n; ++i) {
    P.push_back(a.control_points()[i]);
    P.push_back(b.control_points()[i]);
    W.push_back(a.weights()[i]);
    W.push_back(b.weights()[i]);
  }
  return {a.degree(), 1, a.knots(), {0, 0, 1, 1}, n, 2, std::move(P), std::move(W)};
}

}  // namespace nf
