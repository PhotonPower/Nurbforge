#include "nurbforge/modeling/fillet.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <string>
#include <utility>
#include <vector>

#include "nurbforge/analysis/closest_point.hpp"
#include "nurbforge/analysis/intersect.hpp"
#include "nurbforge/analysis/surface_intersect.hpp"
#include "nurbforge/core/error.hpp"
#include "nurbforge/fitting/fit.hpp"
#include "nurbforge/geometry/primitives.hpp"
#include "nurbforge/modeling/loft.hpp"
#include "nurbforge/modeling/offset.hpp"

namespace nf {

namespace {

std::string num(double x) { return std::to_string(x); }

double sign_of(double x) { return x > 0.0 ? 1.0 : (x < 0.0 ? -1.0 : 0.0); }

// ------------------------------------------------------------------ curves

// Point, unit tangent and signed curvature (relative to the plane normal n) of a planar curve.
struct PlanarFrame {
  Vec3 p, d1, t;
  double kappa{0.0};
};

PlanarFrame planar_frame(const NurbsCurve& c, double s, const Vec3& n) {
  const auto D = c.derivatives(s, 2);
  const double len = norm(D[1]);
  if (len <= tol::kLength)
    throw GeometryError("fillet_curves: the curve tangent vanishes at t = " + num(s));
  return {D[0], D[1], D[1] / len, dot(n, cross(D[1], D[2])) / (len * len * len)};
}

// Offset point c(s) + d * normalize(n x c'(s)) and its derivative c'(s) (1 - d kappa).
std::pair<Vec3, Vec3> offset_point(const NurbsCurve& c, double s, double d, const Vec3& n) {
  const PlanarFrame f = planar_frame(c, s, n);
  return {f.p + d * cross(n, f.t), f.d1 * (1.0 - d * f.kappa)};
}

// Gauss-Newton on a(s) + da Na(s) = b(t) + db Nb(t) (box-constrained to the domains).
void refine_center(const NurbsCurve& a, const NurbsCurve& b, double da, double db, const Vec3& n,
                   double& s, double& t) {
  const double scale = std::max(1.0, std::max(std::abs(da), std::abs(db)));
  for (int it = 0; it < 50; ++it) {
    const auto [pa, ja] = offset_point(a, s, da, n);
    const auto [pb, jb] = offset_point(b, t, db, n);
    const Vec3 F = pa - pb;
    // J = [ja, -jb]; normal equations J^T J x = -J^T F
    const double m11 = dot(ja, ja), m12 = -dot(ja, jb), m22 = dot(jb, jb);
    const double r1 = -dot(ja, F), r2 = dot(jb, F);
    const double det = m11 * m22 - m12 * m12;
    if (!(std::abs(det) > 1e-300)) break;
    const double ds = (r1 * m22 - m12 * r2) / det, dt = (m11 * r2 - m12 * r1) / det;
    s = std::clamp(s + ds, a.domain_start(), a.domain_end());
    t = std::clamp(t + dt, b.domain_start(), b.domain_end());
    if (norm(F) <= 1e-14 * scale && std::abs(ds) + std::abs(dt) <= tol::kParametric) break;
    if (std::abs(ds) * norm(ja) + std::abs(dt) * norm(jb) <= 1e-15 * scale) break;
  }
}

// ------------------------------------------------------------------ surfaces

struct Section {
  double t{0.0};
  Vec3 center, p1, p2;
  double u1{0.0}, v1{0.0}, u2{0.0}, v2{0.0};
  Vec3 x, y;  // unit directions to contact 1 and in the arc plane towards contact 2
  double theta{0.0};
};

class RollingBall {
 public:
  RollingBall(const NurbsSurface& s1, const NurbsSurface& s2, double r,
              const SurfaceFilletOptions& opt, const IntersectionCurve& spine0)
      : s1_(s1), s2_(s2), r_(r), opt_(opt), spine0_(spine0) {}

  // Exact ball centre in the plane through spine0(t) normal to the spine tangent.
  [[nodiscard]] Section solve(double t) const {
    const auto D = spine0_.curve.derivatives(t, 1);
    const Vec3 Q = D[0];
    const double dl = norm(D[1]);
    if (dl <= tol::kLength) throw GeometryError("fillet_surfaces: degenerate spine at t = " + num(t));
    const Vec3 T = D[1] / dl;
    const Vec3 uv1 = spine0_.pcurve1.evaluate(t), uv2 = spine0_.pcurve2.evaluate(t);
    Section sec;
    sec.t = t;
    sec.u1 = uv1.x;
    sec.v1 = uv1.y;
    sec.u2 = uv2.x;
    sec.v2 = uv2.y;
    Vec3 C = Q;
    const double scale = std::max(1.0, r_);
    double res = std::numeric_limits<double>::infinity();
    Vec3 g1, g2;
    const Tolerance ptol{1e-13 * scale, tol::kParametric, tol::kAngular};
    for (int it = 0; it < 40; ++it) {
      const SurfacePoint f1 = closest_point_local(s1_, C, sec.u1, sec.v1, ptol);
      const SurfacePoint f2 = closest_point_local(s2_, C, sec.u2, sec.v2, ptol);
      sec.u1 = f1.u;
      sec.v1 = f1.v;
      sec.u2 = f2.u;
      sec.v2 = f2.v;
      sec.p1 = f1.point;
      sec.p2 = f2.point;
      if (f1.distance <= 1e-9 * scale || f2.distance <= 1e-9 * scale)
        throw GeometryError("fillet_surfaces: the ball centre lies on a surface at t = " + num(t));
      g1 = (C - f1.point) / f1.distance;
      g2 = (C - f2.point) / f2.distance;
      const double e1 = f1.distance - r_, e2 = f2.distance - r_, e3 = dot(C - Q, T);
      res = std::max(std::abs(e1), std::abs(e2));
      // rows g1, g2, T: inverse has columns (g2 x T, T x g1, g1 x g2) / det
      const double det = dot(g1, cross(g2, T));
      if (std::abs(det) <= 1e-12)
        throw GeometryError("fillet_surfaces: the surfaces are tangent at the ball position (t = " +
                            num(t) + ")");
      const Vec3 delta = (-e1 * cross(g2, T) - e2 * cross(T, g1) - e3 * cross(g1, g2)) / det;
      if (res <= 1e-12 * scale && std::abs(e3) <= 1e-12 * scale) break;
      C += delta;
      if (norm(delta) <= 1e-15 * scale) break;
    }
    // Near the end of the spine the exact contact may lie (by the offset fitting error) just
    // outside a parameter domain; the clamped foot point then leaves a small residual.
    if (!(res <= 0.25 * opt_.tolerance))
      throw GeometryError("fillet_surfaces: cannot place the ball at t = " + num(t) +
                          " (residual " + num(res) + ")");
    if (dot(g1, s1_.normal(sec.u1, sec.v1)) * opt_.side1 <= 0.0 ||
        dot(g2, s2_.normal(sec.u2, sec.v2)) * opt_.side2 <= 0.0)
      throw GeometryError("fillet_surfaces: the ball is on the wrong side of a surface at t = " +
                          num(t));
    sec.center = C;
    sec.x = -g1;
    const Vec3 w = -g2;
    sec.theta = std::acos(std::clamp(dot(sec.x, w), -1.0, 1.0));
    if (sec.theta <= 1e-6 || sec.theta >= std::numbers::pi - 1e-6)
      throw GeometryError("fillet_surfaces: degenerate fillet section at t = " + num(t) +
                          " (contact points coincide or are antipodal)");
    sec.y = normalized(w - dot(w, sec.x) * sec.x);
    return sec;
  }

  const NurbsSurface& s1_;
  const NurbsSurface& s2_;
  double r_;
  const SurfaceFilletOptions& opt_;
  const IntersectionCurve& spine0_;
};

// Circular arc of the section with a fixed number of segments (A7.1) on [0, 1], so that all
// sections share one knot vector and u corresponds to the same angle fraction everywhere.
NurbsCurve arc_section(const Section& s, double r, std::size_t nseg) {
  const double dth = s.theta / static_cast<double>(nseg);
  const double w1 = std::cos(dth / 2.0);
  std::vector<Vec3> P{s.center + r * s.x};
  std::vector<double> W{1.0};
  for (std::size_t i = 1; i <= nseg; ++i) {
    const double am = dth * (static_cast<double>(i) - 0.5), ae = dth * static_cast<double>(i);
    P.push_back(s.center + (r / w1) * (std::cos(am) * s.x + std::sin(am) * s.y));
    W.push_back(w1);
    P.push_back(s.center + r * (std::cos(ae) * s.x + std::sin(ae) * s.y));
    W.push_back(1.0);
  }
  return {2, detail::arc_knots(nseg), std::move(P), std::move(W)};
}

double polyline_length(const IntersectionCurve& c) {
  double l = 0.0;
  for (std::size_t i = 1; i < c.points.size(); ++i)
    l += distance(c.points[i - 1].point, c.points[i].point);
  return l;
}

}  // namespace

CurveFillet fillet_curves(const NurbsCurve& a, const NurbsCurve& b, double radius,
                          const Vec3& plane_normal, const CurveFilletOptions& opt) {
  if (!(radius > 0.0)) throw GeometryError("fillet_curves: radius must be positive");
  if (opt.side_a < -1 || opt.side_a > 1 || opt.side_b < -1 || opt.side_b > 1)
    throw GeometryError("fillet_curves: side_a and side_b must be -1, 0 or +1");
  const double nn = norm(plane_normal);
  if (nn <= tol::kLength) throw GeometryError("fillet_curves: plane normal must be non-zero");
  const Vec3 n = plane_normal / nn;

  // corner
  const Vec3 ref = opt.corner ? *opt.corner
                              : 0.5 * (a.evaluate(a.domain_end()) + b.evaluate(b.domain_start()));
  const auto hits = intersect(a, b, Tolerance{1e-9, tol::kParametric, tol::kAngular});
  if (hits.empty()) throw GeometryError("fillet_curves: the curves do not intersect");
  const auto* corner = &hits.front();
  for (const auto& h : hits)
    if (distance(h.point, ref) < distance(corner->point, ref)) corner = &h;
  const double sx = corner->s, tx = corner->t;
  const Vec3 X = corner->point;
  const Vec3 ta = planar_frame(a, sx, n).t, tb = planar_frame(b, tx, n).t;
  if (norm(cross(ta, tb)) <= 1e-9)
    throw GeometryError("fillet_curves: the curves are tangent at the corner");

  // sides (automatic: the sector of the kept parts, which are the longer sides of X)
  const bool keep_a_before = sx - a.domain_start() >= a.domain_end() - sx;
  const bool keep_b_before = tx - b.domain_start() >= b.domain_end() - tx;
  const Vec3 dir_a = keep_a_before ? -ta : ta, dir_b = keep_b_before ? -tb : tb;
  double side_a = opt.side_a, side_b = opt.side_b;
  if (side_a == 0.0) side_a = sign_of(dot(cross(n, ta), dir_b));
  if (side_b == 0.0) side_b = sign_of(dot(cross(n, tb), dir_a));
  if (side_a == 0.0 || side_b == 0.0)
    throw GeometryError("fillet_curves: cannot determine the fillet side");
  const double da = side_a * radius, db = side_b * radius;

  // start value: intersection of the offset curves nearest to the corner
  const NurbsCurve oa = offset_curve(a, da, n, opt.tolerance);
  const NurbsCurve ob = offset_curve(b, db, n, opt.tolerance);
  const auto centers = intersect(oa, ob, Tolerance{std::max(opt.tolerance, 1e-9), tol::kParametric,
                                                   tol::kAngular});
  if (centers.empty())
    throw GeometryError("fillet_curves: the offset curves do not intersect (radius " + num(radius) +
                        " too large for the curves or the chosen sides)");
  const auto* c0 = &centers.front();
  for (const auto& h : centers)
    if (distance(h.point, X) < distance(c0->point, X)) c0 = &h;
  double s = c0->s, t = c0->t;
  refine_center(a, b, da, db, n, s, t);
  const Vec3 ca = offset_point(a, s, da, n).first, cb = offset_point(b, t, db, n).first;
  const double gap = distance(ca, cb);
  if (!(gap <= 1e-9 * std::max(1.0, radius)))
    throw GeometryError("fillet_curves: no fillet centre found (gap " + num(gap) + ")");
  const Vec3 C = 0.5 * (ca + cb), pa = a.evaluate(s), pb = b.evaluate(t);

  // kept parts: the side of the tangent point away from the corner; the arc starts at pa in
  // the direction of travel of the kept part of a (towards the corner)
  const bool a_before = s < sx, b_after = t > tx;
  const Vec3 xa = normalized(pa - C);
  Vec3 ya = planar_frame(a, s, n).t * (a_before ? 1.0 : -1.0);
  ya = normalized(ya - dot(ya, xa) * xa);
  double theta = std::atan2(dot(pb - C, ya), dot(pb - C, xa));
  if (theta <= 0.0) theta += 2.0 * std::numbers::pi;
  CurveFillet out{make_arc(C, xa, ya, radius, 0.0, theta), C, s, t, pa, pb, std::nullopt,
                  std::nullopt};

  if (opt.trim) {
    const double eps = tol::kParametric;
    if (a_before ? s - a.domain_start() <= eps : a.domain_end() - s <= eps)
      throw GeometryError("fillet_curves: the fillet consumes curve a (radius too large)");
    if (b_after ? b.domain_end() - t <= eps : t - b.domain_start() <= eps)
      throw GeometryError("fillet_curves: the fillet consumes curve b (radius too large)");
    out.trimmed_a = a_before ? a.sub_curve(a.domain_start(), s) : a.sub_curve(s, a.domain_end()).reversed();
    out.trimmed_b = b_after ? b.sub_curve(t, b.domain_end()) : b.sub_curve(b.domain_start(), t).reversed();
  }
  return out;
}

SurfaceFillet fillet_surfaces(const NurbsSurface& s1, const NurbsSurface& s2, double radius,
                              const SurfaceFilletOptions& opt) {
  if (!(radius > 0.0)) throw GeometryError("fillet_surfaces: radius must be positive");
  if ((opt.side1 != 1 && opt.side1 != -1) || (opt.side2 != 1 && opt.side2 != -1))
    throw GeometryError("fillet_surfaces: side1 and side2 must be +1 or -1");
  if (!(opt.tolerance > 0.0)) throw GeometryError("fillet_surfaces: tolerance must be positive");
  const double otol = opt.offset_tolerance > 0.0 ? opt.offset_tolerance : 0.1 * opt.tolerance;

  // start values: intersection of the offset surfaces
  const NurbsSurface o1 = offset_surface(s1, opt.side1 * radius, otol);
  const NurbsSurface o2 = offset_surface(s2, opt.side2 * radius, otol);
  const SurfaceIntersection cut = intersect(o1, o2);
  if (cut.curves.empty())
    throw GeometryError("fillet_surfaces: the offset surfaces do not intersect (no room for a "
                        "ball of radius " + num(radius) + " on the chosen sides)");
  const IntersectionCurve* spine0 = &cut.curves.front();
  for (const auto& c : cut.curves)
    if (polyline_length(c) > polyline_length(*spine0)) spine0 = &c;

  const RollingBall ball(s1, s2, radius, opt, *spine0);
  const double t0 = spine0->curve.domain_start(), t1 = spine0->curve.domain_end();

  std::vector<Section> secs;
  constexpr std::size_t kInitial = 9;
  for (std::size_t k = 0; k < kInitial; ++k)
    secs.push_back(ball.solve(k + 1 == kInitial ? t1
                                                : t0 + (t1 - t0) * static_cast<double>(k) /
                                                           static_cast<double>(kInitial - 1)));
  std::size_t nseg = 1;
  for (const auto& sc : secs) nseg = std::max(nseg, detail::arc_segment_count(sc.theta));

  constexpr int kChecksU = 16;
  constexpr double kFractions[] = {0.25, 0.5, 0.75};
  while (true) {
    const std::size_t K = secs.size();
    std::vector<double> ts(K);
    std::vector<NurbsCurve> arcs;
    std::vector<Vec3> centers(K), uv1(K), uv2(K);
    for (std::size_t k = 0; k < K; ++k) {
      ts[k] = secs[k].t;
      arcs.push_back(arc_section(secs[k], radius, nseg));
      centers[k] = secs[k].center;
      uv1[k] = {secs[k].u1, secs[k].v1, 0.0};
      uv2[k] = {secs[k].u2, secs[k].v2, 0.0};
    }
    NurbsSurface surface = loft(arcs, 3, ts);
    NurbsCurve spine = interpolate_curve(centers, 3, ts);
    NurbsCurve pc1 = interpolate_curve(uv1, 3, ts);
    NurbsCurve pc2 = interpolate_curve(uv2, 3, ts);

    std::vector<Section> mids(K - 1);
    std::vector<bool> bad(K - 1, false);
    double worst = 0.0;
    for (std::size_t k = 0; k + 1 < K; ++k) {
      for (double f : kFractions) {
        const double t = ts[k] + f * (ts[k + 1] - ts[k]);
        const Section sc = ball.solve(t);
        const NurbsCurve ideal = arc_section(sc, radius, nseg);
        const Vec3 q1 = pc1.evaluate(t), q2 = pc2.evaluate(t);
        double err = std::max({distance(spine.evaluate(t), sc.center),
                               distance(s1.evaluate(q1.x, q1.y), sc.p1),
                               distance(s2.evaluate(q2.x, q2.y), sc.p2)});
        for (int i = 0; i <= kChecksU; ++i) {
          const double u = i / static_cast<double>(kChecksU);
          err = std::max(err, distance(surface.evaluate(u, t), ideal.evaluate(u)));
        }
        if (f == 0.5) mids[k] = sc;
        worst = std::max(worst, err);
        if (err > opt.tolerance) bad[k] = true;
      }
    }
    const auto nbad = static_cast<std::size_t>(std::count(bad.begin(), bad.end(), true));
    if (nbad == 0) {
      NurbsCurve contact1 = surface.iso_curve_u(0.0);
      NurbsCurve contact2 = surface.iso_curve_u(1.0);
      return {std::move(surface), std::move(spine), std::move(contact1), std::move(contact2),
              std::move(pc1), std::move(pc2)};
    }
    if (K + nbad > opt.max_sections)
      throw GeometryError("fillet_surfaces: tolerance " + num(opt.tolerance) + " not reached with " +
                          std::to_string(K) + " sections (deviation " + num(worst) + ")");
    std::vector<Section> next;
    for (std::size_t k = 0; k < K; ++k) {
      next.push_back(secs[k]);
      if (k + 1 < K && bad[k]) next.push_back(mids[k]);
    }
    secs = std::move(next);
  }
}

}  // namespace nf
