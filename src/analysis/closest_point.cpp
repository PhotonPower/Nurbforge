#include "nurbforge/analysis/closest_point.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace nf {

namespace {

constexpr int kMaxNewton = 60;
constexpr int kMaxDepthCurve = 60;
constexpr int kMaxDepthSurface = 40;
constexpr double kFlatCos = 0.9;  // control polygon edges within ~25° of the main direction

// Backtracking accepts steps that increase |r|^2 only within rounding noise. Near the optimum a
// position error delta changes |r|^2 by ~delta^2, so a strict test would stall at ~sqrt(eps).
constexpr double kNoise = 1.0 + 1e-12;

double dist2(const Vec3& a, const Vec3& b) { return dot(a - b, a - b); }

// Distance from P to the axis-aligned bounding box of the points (lower bound of the distance to
// the curve/surface, by the convex hull property for positive weights).
double box_distance(const std::vector<Vec3>& pts, const Vec3& P) {
  Vec3 lo = pts.front(), hi = pts.front();
  for (const auto& q : pts) {
    lo = {std::min(lo.x, q.x), std::min(lo.y, q.y), std::min(lo.z, q.z)};
    hi = {std::max(hi.x, q.x), std::max(hi.y, q.y), std::max(hi.z, q.z)};
  }
  const Vec3 d{std::max({lo.x - P.x, 0.0, P.x - hi.x}), std::max({lo.y - P.y, 0.0, P.y - hi.y}),
               std::max({lo.z - P.z, 0.0, P.z - hi.z})};
  return norm(d);
}

double box_diagonal(const std::vector<Vec3>& pts) {
  Vec3 lo = pts.front(), hi = pts.front();
  for (const auto& q : pts) {
    lo = {std::min(lo.x, q.x), std::min(lo.y, q.y), std::min(lo.z, q.z)};
    hi = {std::max(hi.x, q.x), std::max(hi.y, q.y), std::max(hi.z, q.z)};
  }
  return distance(lo, hi);
}

// Every edge e satisfies e . dir >= kFlatCos |e| |dir| (zero-length edges allowed).
bool edges_aligned(const std::vector<Vec3>& edges, const Vec3& dir) {
  const double ld = norm(dir);
  if (ld <= 0.0) return false;
  for (const auto& e : edges)
    if (dot(e, dir) < kFlatCos * norm(e) * ld) return false;
  return true;
}

// ------------------------------------------------------------------ curves

CurvePoint newton_curve(const NurbsCurve& c, const Vec3& P, double u, double a, double b,
                        const Tolerance& tol) {
  for (int it = 0; it < kMaxNewton; ++it) {
    const auto D = c.derivatives(u, 2);
    const Vec3 r = D[0] - P;
    const double nr = norm(r), nd = norm(D[1]);
    if (nr <= tol.length || nd == 0.0) break;
    const double f = dot(D[1], r);
    if (f == 0.0) break;  // exactly orthogonal
    const double fp = dot(D[2], r) + nd * nd;
    const double step = fp > 0.0 ? -f / fp : -f / (nd * nd);
    // backtracking: never increase the distance
    double t = 1.0, un = u;
    const double d0 = nr * nr;
    while (true) {
      un = std::clamp(u + t * step, a, b);
      if (dist2(c.evaluate(un), P) <= d0 * kNoise || t < 1e-8) break;
      t *= 0.5;
    }
    const bool done = std::abs(un - u) * nd <= tol.length;
    u = un;
    if (done) break;
  }
  const Vec3 q = c.evaluate(u);
  return {u, q, distance(q, P)};
}

struct CurveSearch {
  const NurbsCurve& curve;
  const Vec3& P;
  const Tolerance& tol;
  CurvePoint best{0.0, {}, std::numeric_limits<double>::max()};

  void consider(const CurvePoint& cp) {
    if (cp.distance < best.distance) best = cp;
  }

  void sample(const NurbsCurve& seg, int n) {
    const double a = seg.domain_start(), b = seg.domain_end();
    for (int i = 0; i <= n; ++i) {
      const double u = a + (b - a) * i / n;
      const Vec3 q = curve.evaluate(u);
      consider({u, q, distance(q, P)});
    }
  }

  bool flat(const NurbsCurve& seg) const {
    const auto& Q = seg.control_points();
    std::vector<Vec3> edges;
    for (std::size_t i = 1; i < Q.size(); ++i) edges.push_back(Q[i] - Q[i - 1]);
    return edges_aligned(edges, Q.back() - Q.front());
  }

  void process(const NurbsCurve& seg, int depth) {
    if (box_distance(seg.control_points(), P) > best.distance + tol.length) return;
    if (depth >= kMaxDepthCurve || flat(seg) || box_diagonal(seg.control_points()) <= tol.length) {
      // start Newton from the best of a few samples inside this segment
      const double a = seg.domain_start(), b = seg.domain_end();
      double u0 = a, d0 = std::numeric_limits<double>::max();
      for (int i = 0; i <= 4; ++i) {
        const double u = a + (b - a) * i / 4.0;
        const double d = dist2(curve.evaluate(u), P);
        if (d < d0) { d0 = d; u0 = u; }
      }
      consider(newton_curve(curve, P, u0, a, b, tol));
      return;
    }
    const double mid = 0.5 * (seg.domain_start() + seg.domain_end());
    const auto [l, r] = seg.split(mid);
    const bool left_first = box_distance(l.control_points(), P) <= box_distance(r.control_points(), P);
    process(left_first ? l : r, depth + 1);
    process(left_first ? r : l, depth + 1);
  }
};

// ------------------------------------------------------------------ surfaces

SurfacePoint newton_surface(const NurbsSurface& s, const Vec3& P, double u, double v,
                            const std::array<double, 2>& ur, const std::array<double, 2>& vr,
                            const Tolerance& tol) {
  for (int it = 0; it < kMaxNewton; ++it) {
    const auto D = s.derivatives(u, v, 2);
    const Vec3 r = D[0][0] - P;
    const Vec3 &Su = D[1][0], &Sv = D[0][1];
    const double nr = norm(r), nu = norm(Su), nv = norm(Sv);
    if (nr <= tol.length) break;
    // gradient of |r|^2 / 2 and its Hessian
    const double f = dot(r, Su), g = dot(r, Sv);
    const double a11 = nu * nu + dot(r, D[2][0]);
    const double a12 = dot(Su, Sv) + dot(r, D[1][1]);
    const double a22 = nv * nv + dot(r, D[0][2]);
    // active set: a coordinate on its bound whose descent direction points outwards stays fixed
    const bool fix_u = (u <= ur[0] && f > 0.0) || (u >= ur[1] && f < 0.0);
    const bool fix_v = (v <= vr[0] && g > 0.0) || (v >= vr[1] && g < 0.0);
    if ((fix_u || f == 0.0) && (fix_v || g == 0.0)) break;  // (projected) stationary point
    const auto newton_1d = [](double grad, double hess, double nt) {
      if (hess > 0.0) return -grad / hess;
      return nt > 0.0 ? -grad / (nt * nt) : 0.0;
    };
    double du = 0.0, dv = 0.0;
    const double det = a11 * a22 - a12 * a12;
    if (fix_u) {
      dv = newton_1d(g, a22, nv);
    } else if (fix_v) {
      du = newton_1d(f, a11, nu);
    } else if (det > 1e-14 * (a11 * a11 + a22 * a22) && a11 > 0.0) {
      du = -(a22 * f - a12 * g) / det;
      dv = -(a11 * g - a12 * f) / det;
    } else {  // indefinite / singular (e.g. pole): scaled gradient step
      du = nu > 0.0 ? -f / (nu * nu) : 0.0;
      dv = nv > 0.0 ? -g / (nv * nv) : 0.0;
    }
    if (du * f + dv * g > 0.0) {  // not a descent direction: fall back to the scaled gradient
      du = fix_u ? 0.0 : (nu > 0.0 ? -f / (nu * nu) : 0.0);
      dv = fix_v ? 0.0 : (nv > 0.0 ? -g / (nv * nv) : 0.0);
    }
    double t = 1.0, un = u, vn = v;
    const double d0 = nr * nr;
    while (true) {
      un = std::clamp(u + t * du, ur[0], ur[1]);
      vn = std::clamp(v + t * dv, vr[0], vr[1]);
      if (dist2(s.evaluate(un, vn), P) <= d0 * kNoise || t < 1e-8) break;
      t *= 0.5;
    }
    const bool done = std::abs(un - u) * nu + std::abs(vn - v) * nv <= tol.length;
    u = un;
    v = vn;
    if (done) break;
  }
  const Vec3 q = s.evaluate(u, v);
  return {u, v, q, distance(q, P)};
}

struct SurfaceSearch {
  const NurbsSurface& surface;
  const Vec3& P;
  const Tolerance& tol;
  SurfacePoint best{0.0, 0.0, {}, std::numeric_limits<double>::max()};

  void consider(const SurfacePoint& sp) {
    if (sp.distance < best.distance) best = sp;
  }

  bool flat(const NurbsSurface& p) const {
    const std::size_t nu = p.num_u(), nv = p.num_v();
    std::vector<Vec3> eu, ev;
    for (std::size_t i = 0; i < nu; ++i)
      for (std::size_t j = 0; j < nv; ++j) {
        if (i + 1 < nu) eu.push_back(p.control_point(i + 1, j) - p.control_point(i, j));
        if (j + 1 < nv) ev.push_back(p.control_point(i, j + 1) - p.control_point(i, j));
      }
    const Vec3 du = (p.control_point(nu - 1, 0) - p.control_point(0, 0)) +
                    (p.control_point(nu - 1, nv - 1) - p.control_point(0, nv - 1));
    const Vec3 dv = (p.control_point(0, nv - 1) - p.control_point(0, 0)) +
                    (p.control_point(nu - 1, nv - 1) - p.control_point(nu - 1, 0));
    return edges_aligned(eu, du) && edges_aligned(ev, dv);
  }

  void process(const NurbsSurface& p, int depth) {
    const auto& Q = p.control_points();
    if (box_distance(Q, P) > best.distance + tol.length) return;
    const auto ur = p.domain_u(), vr = p.domain_v();
    if (depth >= kMaxDepthSurface || flat(p) || box_diagonal(Q) <= tol.length) {
      // Newton from the three best of 5x5 samples
      struct Start { double d, u, v; };
      std::vector<Start> starts;
      for (int i = 0; i <= 4; ++i)
        for (int j = 0; j <= 4; ++j) {
          const double u = ur[0] + (ur[1] - ur[0]) * i / 4.0;
          const double v = vr[0] + (vr[1] - vr[0]) * j / 4.0;
          starts.push_back({dist2(surface.evaluate(u, v), P), u, v});
        }
      std::partial_sort(starts.begin(), starts.begin() + 3, starts.end(),
                        [](const Start& x, const Start& y) { return x.d < y.d; });
      for (int k = 0; k < 3; ++k)
        consider(newton_surface(surface, P, starts[static_cast<std::size_t>(k)].u,
                                starts[static_cast<std::size_t>(k)].v, ur, vr, tol));
      return;
    }
    // split the direction with the larger control-net extent
    double ext_u = 0.0, ext_v = 0.0;
    for (std::size_t j = 0; j < p.num_v(); ++j)
      ext_u = std::max(ext_u, distance(p.control_point(0, j), p.control_point(p.num_u() - 1, j)));
    for (std::size_t i = 0; i < p.num_u(); ++i)
      ext_v = std::max(ext_v, distance(p.control_point(i, 0), p.control_point(i, p.num_v() - 1)));
    const auto [a, b] = ext_u >= ext_v ? p.split_u(0.5 * (ur[0] + ur[1]))
                                       : p.split_v(0.5 * (vr[0] + vr[1]));
    const bool a_first = box_distance(a.control_points(), P) <= box_distance(b.control_points(), P);
    process(a_first ? a : b, depth + 1);
    process(a_first ? b : a, depth + 1);
  }
};

}  // namespace

CurvePoint closest_point(const NurbsCurve& c, const Vec3& P, const Tolerance& tol) {
  CurveSearch search{c, P, tol};
  const auto segments = c.decompose_bezier();
  for (const auto& seg : segments) search.sample(seg, 2);  // initial upper bound
  for (const auto& seg : segments) search.process(seg, 0);
  return search.best;
}

CurvePoint closest_point_local(const NurbsCurve& c, const Vec3& P, double u0, const Tolerance& tol) {
  const double a = c.domain_start(), b = c.domain_end();
  return newton_curve(c, P, std::clamp(u0, a, b), a, b, tol);
}

SurfacePoint closest_point(const NurbsSurface& s, const Vec3& P, const Tolerance& tol) {
  SurfaceSearch search{s, P, tol};
  const auto patches = s.decompose_bezier();
  for (const auto& p : patches) {  // initial upper bound: patch corners and centre
    const auto ur = p.domain_u(), vr = p.domain_v();
    for (double u : {ur[0], 0.5 * (ur[0] + ur[1]), ur[1]})
      for (double v : {vr[0], 0.5 * (vr[0] + vr[1]), vr[1]}) {
        const Vec3 q = s.evaluate(u, v);
        search.consider({u, v, q, distance(q, P)});
      }
  }
  for (const auto& p : patches) search.process(p, 0);
  return search.best;
}

SurfacePoint closest_point_local(const NurbsSurface& s, const Vec3& P, double u0, double v0,
                                 const Tolerance& tol) {
  const auto ur = s.domain_u(), vr = s.domain_v();
  return newton_surface(s, P, std::clamp(u0, ur[0], ur[1]), std::clamp(v0, vr[0], vr[1]), ur, vr, tol);
}

}  // namespace nf
