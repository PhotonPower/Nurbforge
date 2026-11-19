#include "nurbforge/modeling/surface_tools.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <string>
#include <vector>

#include "nurbforge/analysis/intersect.hpp"
#include "nurbforge/core/error.hpp"
#include "nurbforge/fitting/fit.hpp"
#include "nurbforge/geometry/transform.hpp"
#include "nurbforge/modeling/construct.hpp"
#include "nurbforge/modeling/extend.hpp"
#include "nurbforge/modeling/gordon.hpp"
#include "nurbforge/modeling/offset.hpp"
#include "nurbforge/nurbs/join.hpp"
#include "nurbforge/nurbs/knot_vector.hpp"

namespace nf {

namespace {

// Clamped uniform knots for n control points on [a, b].
KnotVector uniform_on(std::size_t degree, std::size_t n, double a, double b) {
  KnotVector k = clamped_uniform_knots(degree, n);
  const double k0 = k.front(), k1 = k.back();
  for (auto& x : k) x = a + (b - a) * (x - k0) / (k1 - k0);
  k.front() = a;
  k.back() = b;
  for (std::size_t i = 0; i <= degree; ++i) {
    k[i] = a;
    k[k.size() - 1 - i] = b;
  }
  return k;
}

std::vector<double> samples(double a, double b, std::size_t n) {
  std::vector<double> t(n);
  for (std::size_t i = 0; i < n; ++i)
    t[i] = i + 1 == n ? b : a + (b - a) * static_cast<double>(i) / static_cast<double>(n - 1);
  return t;
}

Vec3 at_start(const NurbsCurve& c) { return c.evaluate(c.domain_start()); }
Vec3 at_end(const NurbsCurve& c) { return c.evaluate(c.domain_end()); }

// The curve split at its corners (tangent discontinuities at knots).
std::vector<NurbsCurve> smooth_pieces(const NurbsCurve& curve) {
  std::vector<NurbsCurve> out;
  NurbsCurve rest = curve.clamped();
  const auto [knots, mult] = unique_knots(rest.knots());
  std::vector<double> cuts;
  for (std::size_t i = 1; i + 1 < knots.size(); ++i) {
    if (static_cast<std::size_t>(mult[i]) < rest.degree()) continue;  // C1 or smoother
    const double t = knots[i], h = 1e-9 * (rest.domain_end() - rest.domain_start());
    const Vec3 a = normalized(rest.derivatives(t - h, 1)[1]);
    const Vec3 b = normalized(rest.derivatives(t + h, 1)[1]);
    if (norm(cross(a, b)) > 1e-7 || dot(a, b) < 0.0) cuts.push_back(t);
  }
  for (const double t : cuts) {
    auto [p, q] = rest.split(t);
    out.push_back(std::move(p));
    rest = std::move(q);
  }
  out.push_back(std::move(rest));
  // a closed curve that is smooth at its start: the first and last piece are one piece
  if (out.size() > 1 &&
      distance(at_start(curve), at_end(curve)) <= 1e-12 * (1.0 + norm(at_start(curve)))) {
    const double a = curve.domain_start(), b = curve.domain_end(), h = 1e-9 * (b - a);
    const Vec3 ta = normalized(curve.derivatives(a + h, 1)[1]);
    const Vec3 tb = normalized(curve.derivatives(b - h, 1)[1]);
    if (norm(cross(ta, tb)) <= 1e-7 && dot(ta, tb) > 0.0) {
      NurbsCurve first = out.front(), last = out.back();
      first = reparametrize(first, last.domain_end(),
                            last.domain_end() + (first.domain_end() - first.domain_start()));
      out.back() = join(last, first, 1e-9);
      out.erase(out.begin());
    }
  }
  return out;
}

bool is_closed(const NurbsCurve& c) {
  double size = 0.0;
  for (const auto& p : c.control_points())
    size = std::max(size, distance(p, c.control_points()[0]));
  return distance(at_start(c), at_end(c)) <= 1e-9 * std::max(1.0, size);
}

// Offsets of the smooth pieces of a planar curve by d, mitred at the corners: neighbouring
// offsets are extended (linearly) or cut to their intersection. Returns pairs (piece, offset
// piece).
std::vector<std::pair<NurbsCurve, NurbsCurve>> mitred_offsets(const NurbsCurve& curve, double d,
                                                              const Vec3& normal, double tol,
                                                              const std::string& who) {
  const std::vector<NurbsCurve> pieces = smooth_pieces(curve);
  // a closed curve with corners has its seam at a corner (smooth seams are merged away)
  const bool closed = is_closed(curve) && pieces.size() > 1;
  std::vector<NurbsCurve> off;
  for (const auto& p : pieces) off.push_back(offset_curve(p, d, normal, tol));
  const std::size_t n = pieces.size();
  std::vector<double> t0(n), t1(n);  // trimmed ranges on the extended offsets
  std::vector<NurbsCurve> ext;
  const double L = 4.0 * std::abs(d) + 1e-9;
  for (const auto& o : off) {
    NurbsCurve e = extend(o, CurveEnd::Start, L, Extension::Linear);
    e = extend(e, CurveEnd::End, L, Extension::Linear);
    ext.push_back(std::move(e));
  }
  for (std::size_t i = 0; i < n; ++i) {
    t0[i] = off[i].domain_start();
    t1[i] = off[i].domain_end();
  }
  const Tolerance itol{1e-10, tol::kParametric, tol::kAngular};
  for (std::size_t i = 0; i + 1 < n + (closed ? 1 : 0); ++i) {
    const std::size_t j = (i + 1) % n;
    const Vec3 corner = 0.5 * (at_end(off[i]) + at_start(off[j]));
    const auto hits = intersect(ext[i], ext[j], itol);
    if (hits.empty())
      throw GeometryError(who + ": the offsets of two pieces do not meet at a corner");
    const auto* best = &hits.front();
    for (const auto& h : hits)
      if (distance(h.point, corner) < distance(best->point, corner)) best = &h;
    t1[i] = best->s;
    t0[j] = best->t;
  }
  std::vector<std::pair<NurbsCurve, NurbsCurve>> out;
  for (std::size_t i = 0; i < n; ++i) {
    if (!(t1[i] > t0[i])) throw GeometryError(who + ": an offset piece vanishes at the corners");
    out.emplace_back(pieces[i], ext[i].sub_curve(t0[i], t1[i]));
  }
  return out;
}

}  // namespace

RebuiltCurve rebuild_curve(const NurbsCurve& curve, std::size_t num_ctrl, std::size_t degree) {
  if (degree < 1 || num_ctrl < degree + 1)
    throw GeometryError("rebuild_curve: need degree >= 1 and num_ctrl >= degree + 1");
  const double a = curve.domain_start(), b = curve.domain_end();
  const std::vector<double> t = samples(a, b, std::max<std::size_t>(10 * num_ctrl, 200));
  std::vector<Vec3> pts;
  for (const double x : t) pts.push_back(curve.evaluate(x));
  NurbsCurve r = approximate_curve(pts, degree, uniform_on(degree, num_ctrl, a, b), t);
  double err = 0.0;
  for (const double x : samples(a, b, 4 * t.size()))
    err = std::max(err, distance(curve.evaluate(x), r.evaluate(x)));
  return {std::move(r), err};
}

RebuiltSurface rebuild_surface(const NurbsSurface& s, std::size_t num_u, std::size_t num_v,
                               std::size_t degree_u, std::size_t degree_v) {
  if (degree_u < 1 || degree_v < 1 || num_u < degree_u + 1 || num_v < degree_v + 1)
    throw GeometryError("rebuild_surface: need degrees >= 1 and num >= degree + 1");
  const auto du = s.domain_u(), dv = s.domain_v();
  const std::vector<double> tu = samples(du[0], du[1], std::max<std::size_t>(6 * num_u, 60));
  const std::vector<double> tv = samples(dv[0], dv[1], std::max<std::size_t>(6 * num_v, 60));
  const KnotVector ku = uniform_on(degree_u, num_u, du[0], du[1]);
  const KnotVector kv = uniform_on(degree_v, num_v, dv[0], dv[1]);
  // rows along v (one per sample u), then columns of their control points along u
  std::vector<NurbsCurve> rows;
  for (const double u : tu) {
    std::vector<Vec3> pts;
    for (const double v : tv) pts.push_back(s.evaluate(u, v));
    rows.push_back(approximate_curve(pts, degree_v, kv, tv));
  }
  std::vector<Vec3> P(num_u * num_v);
  for (std::size_t j = 0; j < num_v; ++j) {
    std::vector<Vec3> col;
    for (const auto& r : rows) col.push_back(r.control_points()[j]);
    const NurbsCurve c = approximate_curve(col, degree_u, ku, tu);
    for (std::size_t i = 0; i < num_u; ++i) P[i * num_v + j] = c.control_points()[i];
  }
  NurbsSurface r(degree_u, degree_v, ku, kv, num_u, num_v, std::move(P));
  double err = 0.0;
  for (const double u : samples(du[0], du[1], 4 * tu.size()))
    for (const double v : samples(dv[0], dv[1], 4 * tv.size()))
      err = std::max(err, distance(s.evaluate(u, v), r.evaluate(u, v)));
  return {std::move(r), err};
}

std::vector<NurbsSurface> extrude_tapered(const NurbsCurve& profile, const Vec3& direction,
                                          double angle, double tolerance) {
  const double h = norm(direction);
  if (h <= tol::kLength) throw GeometryError("extrude_tapered: direction must be non-zero");
  if (!(std::abs(angle) < 0.5 * std::numbers::pi))
    throw GeometryError("extrude_tapered: |angle| must be below 90 degrees");
  const Vec3 n = direction / h;
  // offset side: outwards = right of a counter-clockwise closed profile (seen from n), left for
  // an open one
  double side = 1.0;
  if (is_closed(profile)) {
    double area = 0.0;
    const std::vector<double> t = samples(profile.domain_start(), profile.domain_end(), 257);
    for (std::size_t i = 0; i + 1 < t.size(); ++i)
      area += dot(n, cross(profile.evaluate(t[i]), profile.evaluate(t[i + 1])));
    side = area > 0.0 ? -1.0 : 1.0;
  }
  const double d = side * h * std::tan(angle);
  std::vector<NurbsSurface> out;
  if (d == 0.0) {
    for (const auto& p : smooth_pieces(profile)) out.push_back(extrude(p, direction));
    return out;
  }
  for (const auto& [piece, top] : mitred_offsets(profile, d, n, tolerance, "extrude_tapered"))
    out.push_back(
        ruled(reparametrize(piece, 0.0, 1.0),
              reparametrize(top, 0.0, 1.0).transformed(Transform::translation(direction))));
  return out;
}

std::vector<NurbsSurface> ribbon(const NurbsCurve& curve, double width, const Vec3& plane_normal,
                                 double tolerance) {
  if (width == 0.0) throw GeometryError("ribbon: width must be non-zero");
  const double nn = norm(plane_normal);
  if (nn <= tol::kLength) throw GeometryError("ribbon: plane normal must be non-zero");
  std::vector<NurbsSurface> out;
  for (const auto& [piece, side] :
       mitred_offsets(curve, width, plane_normal / nn, tolerance, "ribbon"))
    out.push_back(ruled(reparametrize(piece, 0.0, 1.0), reparametrize(side, 0.0, 1.0)));
  return out;
}

NurbsSurface extrude_to_point(const NurbsCurve& curve, const Vec3& apex) {
  // the apex as a degenerate curve with the curve's degree, knots and weights: the rulings run
  // from c(u) to the apex in homogeneous coordinates, which is exact also for rational curves
  const std::vector<Vec3> P(curve.num_control_points(), apex);
  return ruled(curve, NurbsCurve(curve.degree(), curve.knots(), P, curve.weights()));
}

NurbsSurface triangle_patch(const NurbsCurve& c0, const NurbsCurve& c1, const NurbsCurve& c2,
                            double tolerance) {
  // A = start of c0, B = end of c0, C = the third corner; d0 joins A and C, d1 joins B and C
  const Vec3 A = at_start(c0), B = at_end(c0);
  const auto touches = [&](const NurbsCurve& c, const Vec3& p) {
    return std::min(distance(at_start(c), p), distance(at_end(c), p));
  };
  const bool c1_at_A = touches(c1, A) + touches(c2, B) < touches(c1, B) + touches(c2, A);
  const NurbsCurve& d0 = c1_at_A ? c1 : c2;
  const NurbsCurve& d1 = c1_at_A ? c2 : c1;
  const Vec3 C = distance(at_start(d0), A) < distance(at_end(d0), A) ? at_end(d0) : at_start(d0);
  const double gap = std::max({touches(d0, A), touches(d1, B), touches(d1, C), touches(d0, C)});
  if (gap > tolerance)
    throw GeometryError("triangle_patch: the curves do not form a closed loop (corner gap " +
                        std::to_string(gap) + ")");
  const NurbsCurve apex(1, {0.0, 0.0, 1.0, 1.0}, {C, C});
  return coons(c0, apex, d0, d1, tolerance);
}

}  // namespace nf
