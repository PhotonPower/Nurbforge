#include "nurbforge/modeling/extend.hpp"

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include "nurbforge/analysis/properties.hpp"
#include "nurbforge/core/error.hpp"
#include "nurbforge/nurbs/join.hpp"

namespace nf {

namespace {

NurbsCurve shift(const NurbsCurve& c, double d) {
  KnotVector U = c.knots();
  for (double& k : U) k += d;
  return {c.degree(), std::move(U), c.control_points(), c.weights()};
}

NurbsSurface shift_u(const NurbsSurface& s, double d) {
  KnotVector U = s.knots_u();
  for (double& k : U) k += d;
  return {s.degree_u(), s.degree_v(), std::move(U), s.knots_v(), s.num_u(), s.num_v(),
          s.control_points(), s.weights()};
}

// Same surface with the u direction reversed: r(u, v) = s(a + b - u, v).
NurbsSurface reverse_u(const NurbsSurface& s) {
  const auto d = s.domain_u();
  const KnotVector& K = s.knots_u();
  KnotVector U(K.size());
  for (std::size_t i = 0; i < K.size(); ++i) U[i] = d[0] + d[1] - K[K.size() - 1 - i];
  std::vector<Vec3> P;
  std::vector<double> W;
  for (std::size_t i = s.num_u(); i-- > 0;)
    for (std::size_t j = 0; j < s.num_v(); ++j) {
      P.push_back(s.control_point(i, j));
      W.push_back(s.weight(i, j));
    }
  return {s.degree_u(), s.degree_v(), std::move(U), s.knots_v(), s.num_u(), s.num_v(),
          std::move(P), std::move(W)};
}

// Clamped curve whose last knot span is a Bezier piece (its start knot raised to multiplicity p).
struct Bezierized {
  NurbsCurve curve;
  double a_last;          // start of the last span
  std::size_t inserted;   // copies of a_last inserted
};

Bezierized bezierize_end(const NurbsCurve& c) {
  const std::size_t p = c.degree();
  const KnotVector& U = c.knots();
  const double b = c.domain_end();
  std::size_t k = c.num_control_points();
  while (U[k] >= b) --k;
  const double a_last = U[k];
  const bool first_span = a_last == c.domain_start();
  const std::size_t mult = knot_multiplicity(U, a_last);
  const std::size_t inserted = first_span || mult >= p ? 0 : p - mult;
  if (inserted == 0) return {c, a_last, 0};
  return {c.refine_knots(std::vector<double>(inserted, a_last)), a_last, inserted};
}

// Natural extension of the last Bezier piece to [a_last, b + delta] (de Casteljau at t > 1, on
// homogeneous points). Knots are not cleaned up.
NurbsCurve natural_raw(const Bezierized& bz, double delta) {
  const NurbsCurve& c = bz.curve;
  const std::size_t p = c.degree();
  auto Pw = c.homogeneous_control_points();
  KnotVector U = c.knots();
  const double b = c.domain_end();
  const double t = (b + delta - bz.a_last) / (b - bz.a_last);
  std::vector<Vec4> tmp(Pw.end() - static_cast<std::ptrdiff_t>(p + 1), Pw.end());
  std::vector<Vec4> L(p + 1);
  L[0] = tmp[0];
  for (std::size_t r = 1; r <= p; ++r) {
    for (std::size_t i = 0; i + r <= p; ++i) tmp[i] = (1.0 - t) * tmp[i] + t * tmp[i + 1];
    L[r] = tmp[0];
  }
  double wmax = 0.0;
  for (const auto& q : Pw) wmax = std::max(wmax, q.w);
  for (const auto& q : L)
    if (!(q.w > 1e-12 * wmax))
      throw GeometryError("extend: natural extension of this rational curve is too long (weights "
                          "would become non-positive); use a shorter length or Extension::Linear");
  std::copy(L.begin(), L.end(), Pw.end() - static_cast<std::ptrdiff_t>(p + 1));
  std::fill(U.end() - static_cast<std::ptrdiff_t>(p + 1), U.end(), b + delta);
  return NurbsCurve::from_homogeneous(p, std::move(U), Pw);
}

// Parameter amount delta whose natural extension has arc length `length` (safeguarded Newton).
double natural_delta(const Bezierized& bz, double length) {
  const double b = bz.curve.domain_end();
  const double speed = norm(bz.curve.derivatives(b, 1)[1]);
  if (speed <= tol::kLength) throw GeometryError("extend: vanishing end tangent");
  double delta = length / speed;
  for (int it = 0; it < 60; ++it) {
    const NurbsCurve e = natural_raw(bz, delta);
    const double g = arc_length(e, b, b + delta) - length;
    if (std::abs(g) <= 1e-12 * length) break;
    const double sp = norm(e.derivatives(b + delta, 1)[1]);
    double next = sp > 0.0 ? delta - g / sp : 0.5 * delta;
    if (!(next > 0.0)) next = 0.5 * delta;
    delta = next;
  }
  return delta;
}

NurbsCurve extend_end(const NurbsCurve& c_in, double length, Extension type) {
  const NurbsCurve c = c_in.clamped();
  const double b = c.domain_end();
  if (type == Extension::Natural) {
    const auto bz = bezierize_end(c);
    const NurbsCurve e = natural_raw(bz, natural_delta(bz, length));
    return bz.inserted > 0 ? e.try_remove_knot(bz.a_last, bz.inserted).first : e;
  }
  const Vec3 D = c.derivatives(b, 1)[1];
  const double speed = norm(D);
  if (speed <= tol::kLength) throw GeometryError("extend: vanishing end tangent");
  const Vec3 P = c.control_points().back();  // = C(b), exactly
  const double delta = length / speed;       // matches the parameter speed: C1
  NurbsCurve line{1, {b, b, b + delta, b + delta}, {P, P + (length / speed) * D}, {}};
  if (c.degree() > 1) line = line.elevate_degree(c.degree() - 1);
  return join(c, line);
}

NurbsSurface surface_from_columns(const NurbsSurface& s, std::size_t pu,
                                  const std::vector<NurbsCurve>& cols) {
  const std::size_t nu = cols.front().num_control_points(), nv = cols.size();
  std::vector<Vec3> P(nu * nv);
  std::vector<double> W(nu * nv);
  for (std::size_t j = 0; j < nv; ++j)
    for (std::size_t i = 0; i < nu; ++i) {
      P[i * nv + j] = cols[j].control_points()[i];
      W[i * nv + j] = cols[j].weights()[i];
    }
  return {pu, s.degree_v(), cols.front().knots(), s.knots_v(), nu, nv, std::move(P), std::move(W)};
}

NurbsSurface extend_umax(const NurbsSurface& s_in, double length, Extension type) {
  const NurbsSurface s = s_in.clamped();
  const std::size_t pu = s.degree_u(), nu = s.num_u(), nv = s.num_v();
  const double b = s.domain_u()[1];
  const auto dv = s.domain_v();
  const double vmid = 0.5 * (dv[0] + dv[1]);

  if (type == Extension::Natural) {
    std::vector<Bezierized> bz;
    for (std::size_t j = 0; j < nv; ++j) {
      std::vector<Vec3> P(nu);
      std::vector<double> W(nu);
      for (std::size_t i = 0; i < nu; ++i) {
        P[i] = s.control_point(i, j);
        W[i] = s.weight(i, j);
      }
      bz.push_back(bezierize_end(NurbsCurve(pu, s.knots_u(), std::move(P), std::move(W))));
    }
    const double delta = natural_delta(bezierize_end(s.iso_curve_v(vmid)), length);
    std::vector<NurbsCurve> cols;
    for (const auto& z : bz) cols.push_back(natural_raw(z, delta));
    const NurbsSurface e = surface_from_columns(s, pu, cols);
    return bz.front().inserted > 0 ? e.try_remove_knot_u(bz.front().a_last, bz.front().inserted).first : e;
  }

  // linear: ruled extension along the boundary tangents, C1. Needs non-rational iso-curves in
  // the extension direction (constant weight along every column); weights may vary across.
  for (std::size_t j = 0; j < nv; ++j)
    for (std::size_t i = 0; i < nu; ++i)
      if (std::abs(s.weight(i, j) - s.weight(0, j)) > tol::kLength * s.weight(0, j))
        throw GeometryError("extend: linear extension needs constant weights along the extension "
                            "direction; use Extension::Natural");
  const double speed = norm(s.derivatives(b, vmid, 1)[1][0]);
  if (speed <= tol::kLength) throw GeometryError("extend: vanishing boundary tangent");
  const double delta = length / speed;
  const double factor = static_cast<double>(pu) / (b - s.knots_u()[nu - 1]);  // C'(b) = factor (P_n - P_n-1)
  std::vector<Vec3> P;
  std::vector<double> Wext;
  for (std::size_t j = 0; j < nv; ++j) {
    P.push_back(s.control_point(nu - 1, j));
    Wext.push_back(s.weight(nu - 1, j));
  }
  for (std::size_t j = 0; j < nv; ++j) {
    P.push_back(s.control_point(nu - 1, j) +
                delta * factor * (s.control_point(nu - 1, j) - s.control_point(nu - 2, j)));
    Wext.push_back(s.weight(nu - 1, j));
  }
  NurbsSurface ext{1, s.degree_v(), {b, b, b + delta, b + delta}, s.knots_v(), 2, nv,
                   std::move(P), std::move(Wext)};
  if (pu > 1) ext = ext.elevate_degree_u(pu - 1);
  return join_u(s, ext);
}

NurbsSurface extend_umin(const NurbsSurface& s, double length, Extension type) {
  const double b = s.domain_u()[1];
  const NurbsSurface e = extend_umax(reverse_u(s), length, type);
  const double delta = e.domain_u()[1] - b;
  return shift_u(reverse_u(e), -delta);
}

}  // namespace

NurbsCurve extend(const NurbsCurve& c, CurveEnd end, double length, Extension type) {
  if (!(length > 0.0)) throw GeometryError("extend: length must be positive");
  if (end == CurveEnd::End) return extend_end(c, length, type);
  const double b = c.domain_end();
  const NurbsCurve e = extend_end(c.reversed(), length, type);
  return shift(e.reversed(), -(e.domain_end() - b));
}

NurbsSurface extend(const NurbsSurface& s, SurfaceSide side, double length, Extension type) {
  if (!(length > 0.0)) throw GeometryError("extend: length must be positive");
  switch (side) {
    case SurfaceSide::UMax: return extend_umax(s, length, type);
    case SurfaceSide::UMin: return extend_umin(s, length, type);
    case SurfaceSide::VMax: return extend_umax(s.transposed(), length, type).transposed();
    case SurfaceSide::VMin: return extend_umin(s.transposed(), length, type).transposed();
  }
  throw GeometryError("extend: unknown side");
}

}  // namespace nf
