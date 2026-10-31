#include "nurbforge/nurbs/join.hpp"

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include "nurbforge/core/error.hpp"

namespace nf {

namespace {

// Knots of `from` (with multiplicity) missing in `into`.
std::vector<double> missing_knots(const KnotVector& from, const KnotVector& into) {
  std::vector<double> X;
  const auto [vals, mults] = unique_knots(from);
  for (std::size_t i = 0; i < vals.size(); ++i)
    for (auto c = knot_multiplicity(into, vals[i]); c < static_cast<std::size_t>(mults[i]); ++c)
      X.push_back(vals[i]);
  return X;
}

// Linear reparametrisation of the v direction to [a, b].
NurbsSurface reparametrize_v(const NurbsSurface& s, double a, double b) {
  const auto d = s.domain_v();
  KnotVector kv = s.knots_v();
  for (double& x : kv) x = a + (x - d[0]) * (b - a) / (d[1] - d[0]);
  // s is clamped: make the end knots exact so rounding cannot break monotonicity
  for (std::size_t i = 0; i <= s.degree_v(); ++i) kv[i] = a;
  for (std::size_t i = s.num_v(); i < kv.size(); ++i) kv[i] = b;
  return {s.degree_u(), s.degree_v(), s.knots_u(), std::move(kv), s.num_u(), s.num_v(),
          s.control_points(), s.weights()};
}

std::string gap_message(double gap, double tolerance) {
  return "join: boundaries do not meet (gap " + std::to_string(gap) + " > tolerance " +
         std::to_string(tolerance) + ")";
}

}  // namespace

NurbsCurve join(const NurbsCurve& a_in, const NurbsCurve& b_in, double tolerance) {
  NurbsCurve a = a_in.clamped(), b = b_in.clamped();
  if (a.degree() < b.degree()) a = a.elevate_degree(b.degree() - a.degree());
  if (b.degree() < a.degree()) b = b.elevate_degree(a.degree() - b.degree());
  const std::size_t p = a.degree();

  const Vec3 pa = a.control_points().back(), pb = b.control_points().front();
  if (distance(pa, pb) > tolerance) throw GeometryError(gap_message(distance(pa, pb), tolerance));

  auto P = a.homogeneous_control_points();
  auto Pb = b.homogeneous_control_points();
  const double scale = P.back().w / Pb.front().w;  // weights are scale-invariant per curve
  P.back() = Vec4::from_point(0.5 * (pa + pb), P.back().w);
  for (std::size_t i = 1; i < Pb.size(); ++i) P.push_back(scale * Pb[i]);

  const double joint = a.domain_end();
  const double shift = joint - b.domain_start();
  KnotVector U(a.knots().begin(), a.knots().end() - 1);  // joint keeps multiplicity p
  for (std::size_t i = p + 1; i < b.knots().size(); ++i) U.push_back(b.knots()[i] + shift);

  const NurbsCurve c = NurbsCurve::from_homogeneous(p, std::move(U), P);
  return c.try_remove_knot(joint, p).first;
}

NurbsSurface join_u(const NurbsSurface& a_in, const NurbsSurface& b_in, double tolerance) {
  NurbsSurface a = a_in.clamped(), b = b_in.clamped();
  // u degree
  if (a.degree_u() < b.degree_u()) a = a.elevate_degree_u(b.degree_u() - a.degree_u());
  if (b.degree_u() < a.degree_u()) b = b.elevate_degree_u(a.degree_u() - b.degree_u());
  // v compatibility: degree, domain, knots
  if (a.degree_v() < b.degree_v()) a = a.elevate_degree_v(b.degree_v() - a.degree_v());
  if (b.degree_v() < a.degree_v()) b = b.elevate_degree_v(a.degree_v() - b.degree_v());
  const auto dv = a.domain_v();
  b = reparametrize_v(b, dv[0], dv[1]);
  const auto Xa = missing_knots(b.knots_v(), a.knots_v());
  const auto Xb = missing_knots(a.knots_v(), b.knots_v());
  a = a.refine_knots_v(Xa);
  b = b.refine_knots_v(Xb);

  const std::size_t pu = a.degree_u(), nv = a.num_v(), nua = a.num_u(), nub = b.num_u();
  // boundary rows: positions and weight ratio
  double gap = 0.0;
  for (std::size_t j = 0; j < nv; ++j)
    gap = std::max(gap, distance(a.control_point(nua - 1, j), b.control_point(0, j)));
  if (gap > tolerance) throw GeometryError(gap_message(gap, tolerance));
  const double scale = a.weight(nua - 1, 0) / b.weight(0, 0);
  for (std::size_t j = 1; j < nv; ++j)
    if (std::abs(a.weight(nua - 1, j) / b.weight(0, j) - scale) > tol::kLength * scale)
      throw GeometryError("join: boundary weights of the two surfaces are not proportional");

  std::vector<Vec3> P;
  std::vector<double> W;
  for (std::size_t i = 0; i < nua; ++i)
    for (std::size_t j = 0; j < nv; ++j) {
      const bool seam = i == nua - 1;
      P.push_back(seam ? 0.5 * (a.control_point(i, j) + b.control_point(0, j)) : a.control_point(i, j));
      W.push_back(a.weight(i, j));
    }
  for (std::size_t i = 1; i < nub; ++i)
    for (std::size_t j = 0; j < nv; ++j) {
      P.push_back(b.control_point(i, j));
      W.push_back(scale * b.weight(i, j));
    }

  const double seam = a.domain_u()[1];
  const double shift = seam - b.domain_u()[0];
  KnotVector U(a.knots_u().begin(), a.knots_u().end() - 1);
  for (std::size_t i = pu + 1; i < b.knots_u().size(); ++i) U.push_back(b.knots_u()[i] + shift);

  const NurbsSurface s{pu, a.degree_v(), std::move(U), a.knots_v(), nua + nub - 1, nv,
                       std::move(P), std::move(W)};
  return s.try_remove_knot_u(seam, pu).first;
}

NurbsSurface join_v(const NurbsSurface& a, const NurbsSurface& b, double tolerance) {
  return join_u(a.transposed(), b.transposed(), tolerance).transposed();
}

}  // namespace nf
