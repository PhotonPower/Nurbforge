#include "nurbforge/analysis/surface_intersect.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numbers>
#include <optional>

#include "nurbforge/analysis/intersect.hpp"
#include "nurbforge/fitting/fit.hpp"
#include "nurbforge/nurbs/join.hpp"

namespace nf {

namespace {

constexpr double pi = std::numbers::pi;
using X4 = std::array<double, 4>;  // (u1, v1, u2, v2)

// ------------------------------------------------------------------ small linear algebra

template <std::size_t N>
bool solve_n(std::array<std::array<double, N>, N> A, std::array<double, N>& b) {
  for (std::size_t k = 0; k < N; ++k) {
    std::size_t p = k;
    for (std::size_t i = k + 1; i < N; ++i)
      if (std::abs(A[i][k]) > std::abs(A[p][k])) p = i;
    if (std::abs(A[p][k]) < 1e-300) return false;
    std::swap(A[k], A[p]);
    std::swap(b[k], b[p]);
    for (std::size_t i = k + 1; i < N; ++i) {
      const double f = A[i][k] / A[k][k];
      for (std::size_t j = k; j < N; ++j) A[i][j] -= f * A[k][j];
      b[i] -= f * b[k];
    }
  }
  for (std::size_t i = N; i-- > 0;) {
    for (std::size_t j = i + 1; j < N; ++j) b[i] -= A[i][j] * b[j];
    b[i] /= A[i][i];
  }
  return true;
}

double comp(const Vec3& v, std::size_t i) { return i == 0 ? v.x : (i == 1 ? v.y : v.z); }

// ------------------------------------------------------------------ geometry helpers

struct Box3 {
  Vec3 lo, hi;
};
Box3 box(const std::vector<Vec3>& pts) {
  Box3 b{pts.front(), pts.front()};
  for (const auto& p : pts) {
    b.lo = {std::min(b.lo.x, p.x), std::min(b.lo.y, p.y), std::min(b.lo.z, p.z)};
    b.hi = {std::max(b.hi.x, p.x), std::max(b.hi.y, p.y), std::max(b.hi.z, p.z)};
  }
  return b;
}
bool overlap(const Box3& a, const Box3& b, double eps) {
  return a.lo.x <= b.hi.x + eps && b.lo.x <= a.hi.x + eps && a.lo.y <= b.hi.y + eps &&
         b.lo.y <= a.hi.y + eps && a.lo.z <= b.hi.z + eps && b.lo.z <= a.hi.z + eps;
}

struct Cone {
  Vec3 axis;
  double half{pi};
};
Cone edge_cone(const std::vector<Vec3>& edges, const Vec3& axis) {
  const double la = norm(axis);
  if (la <= 0.0) return {};
  Cone c{axis / la, 0.0};
  for (const auto& e : edges) {
    const double le = norm(e);
    if (le > 0.0) c.half = std::max(c.half, std::acos(std::clamp(dot(e, c.axis) / le, -1.0, 1.0)));
  }
  return c;
}
// Conservative estimate of the normal cone of a Bezier patch from its u- and v-edge cones.
Cone normal_cone(const NurbsSurface& s) {
  const std::size_t nu = s.num_u(), nv = s.num_v();
  std::vector<Vec3> eu, ev;
  for (std::size_t i = 0; i < nu; ++i)
    for (std::size_t j = 0; j < nv; ++j) {
      if (i + 1 < nu) eu.push_back(s.control_point(i + 1, j) - s.control_point(i, j));
      if (j + 1 < nv) ev.push_back(s.control_point(i, j + 1) - s.control_point(i, j));
    }
  const Vec3 du = (s.control_point(nu - 1, 0) - s.control_point(0, 0)) +
                  (s.control_point(nu - 1, nv - 1) - s.control_point(0, nv - 1));
  const Vec3 dv = (s.control_point(0, nv - 1) - s.control_point(0, 0)) +
                  (s.control_point(nu - 1, nv - 1) - s.control_point(nu - 1, 0));
  const Cone cu = edge_cone(eu, du), cv = edge_cone(ev, dv);
  const Vec3 n = cross(du, dv);
  if (norm(n) == 0.0 || cu.half >= pi / 2 || cv.half >= pi / 2) return {};
  return {n / norm(n), (cu.half + cv.half) / std::max(norm(cross(cu.axis, cv.axis)), 0.1)};
}

double segment_distance(const Vec3& p, const Vec3& a, const Vec3& b) {
  const Vec3 ab = b - a;
  const double l2 = dot(ab, ab);
  const double t = l2 > 0.0 ? std::clamp(dot(p - a, ab) / l2, 0.0, 1.0) : 0.0;
  return distance(p, a + t * ab);
}

// ------------------------------------------------------------------ the solver

class Ssi {
 public:
  Ssi(const NurbsSurface& a, const NurbsSurface& b, const SurfaceIntersectionOptions& opt)
      : a_(a), b_(b), tol_(opt.tol) {
    std::vector<Vec3> all = a.control_points();
    all.insert(all.end(), b.control_points().begin(), b.control_points().end());
    const Box3 bx = box(all);
    size_ = distance(bx.lo, bx.hi);
    chord_ = opt.chord_tolerance > 0.0 ? opt.chord_tolerance : 1e-4 * size_;
    fit_ = opt.fit_tolerance > 0.0 ? opt.fit_tolerance : 1e-7 * size_;
    hmax_ = size_ / 20.0;
    hmin_ = std::max(1e-9 * size_, 10.0 * tol_.length);
    const auto ua = a.domain_u(), va = a.domain_v(), ub = b.domain_u(), vb = b.domain_v();
    lo_ = {ua[0], va[0], ub[0], vb[0]};
    hi_ = {ua[1], va[1], ub[1], vb[1]};
    const std::array<const KnotVector*, 4> kv{&a.knots_u(), &a.knots_v(), &b.knots_u(), &b.knots_v()};
    for (std::size_t k = 0; k < 4; ++k)
      for (double kn : unique_knots(*kv[k]).first)
        if (kn > lo_[k] && kn < hi_[k]) knots_[k].push_back(kn);
  }

  SurfaceIntersection run() {
    std::vector<X4> seeds;
    const auto pa = a_.decompose_bezier(), pb = b_.decompose_bezier();
    for (const auto& p : pa)
      for (const auto& q : pb) collect_seeds(p, q, 0, seeds);
    // seeds on a domain boundary first: open curves are then traced from their ends
    std::stable_sort(seeds.begin(), seeds.end(), [&](const X4& x, const X4& y) {
      return on_boundary(x) && !on_boundary(y);
    });
    for (const auto& s : seeds) trace_from(s);
    return std::move(result_);
  }

 private:
  // ---------------------------------------------------------------- evaluation
  struct Eval {
    Vec3 p1, p2, su1, sv1, su2, sv2;
  };
  Eval eval(const X4& x) const {
    const auto d1 = a_.first_derivatives(x[0], x[1]);
    const auto d2 = b_.first_derivatives(x[2], x[3]);
    return {d1[0], d2[0], d1[1], d1[2], d2[1], d2[2]};
  }
  X4 clamp(X4 x) const {
    for (std::size_t i = 0; i < 4; ++i) x[i] = std::clamp(x[i], lo_[i], hi_[i]);
    return x;
  }
  bool on_boundary(const X4& x) const {
    for (std::size_t i = 0; i < 4; ++i) {
      const double eps = 1e-12 * (hi_[i] - lo_[i]);
      if (x[i] <= lo_[i] + eps || x[i] >= hi_[i] - eps) return true;
    }
    return false;
  }
  // Unit tangent n1 x n2 (limit normals at poles); nullopt if the surfaces are tangent.
  std::optional<Vec3> tangent(const X4& x) const {
    const Vec3 n1 = a_.normal(x[0], x[1]), n2 = b_.normal(x[2], x[3]);
    const Vec3 t = cross(n1, n2);
    const double l = norm(t);
    if (l <= 1e-7) return std::nullopt;
    return t / l;
  }
  // Parameter velocity (du, dv) with Su du + Sv dv = t in the least-squares sense.
  static std::array<double, 2> param_dir(const Vec3& su, const Vec3& sv, const Vec3& t) {
    const double E = dot(su, su), F = dot(su, sv), G = dot(sv, sv), det = E * G - F * F;
    if (std::abs(det) < 1e-300) return {0.0, 0.0};
    const double a = dot(su, t), b = dot(sv, t);
    return {(G * a - F * b) / det, (E * b - F * a) / det};
  }

  // ---------------------------------------------------------------- Newton solvers
  // Minimum-norm Gauss-Newton on S1 - S2 = 0 (3 equations, 4 unknowns).
  std::optional<X4> refine(X4 x) const {
    for (int it = 0; it < 50; ++it) {
      const Eval e = eval(x);
      const Vec3 r = e.p1 - e.p2;
      if (norm(r) <= 1e-3 * tol_.length) return x;
      const std::array<Vec3, 4> J{e.su1, e.sv1, -e.su2, -e.sv2};
      std::array<std::array<double, 3>, 3> JJ{};
      for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j)
          for (const auto& c : J) JJ[i][j] += comp(c, i) * comp(c, j);
      double tr = JJ[0][0] + JJ[1][1] + JJ[2][2];
      for (std::size_t i = 0; i < 3; ++i) JJ[i][i] += 1e-13 * tr;
      std::array<double, 3> y{-r.x, -r.y, -r.z};
      if (!solve_n<3>(JJ, y)) return std::nullopt;
      const Vec3 yv{y[0], y[1], y[2]};
      X4 xn = x;
      double t = 1.0;
      const double r2 = dot(r, r);
      while (true) {
        for (std::size_t k = 0; k < 4; ++k) xn[k] = x[k] + t * dot(J[k], yv);
        xn = clamp(xn);
        const Eval en = eval(xn);
        if (dot(en.p1 - en.p2, en.p1 - en.p2) <= r2 * (1.0 + 1e-12) || t < 1e-8) break;
        t *= 0.5;
      }
      double step = 0.0;
      for (std::size_t k = 0; k < 4; ++k) step += std::abs(xn[k] - x[k]) * norm(J[k]);
      x = xn;
      if (step <= 1e-15 * size_) break;
    }
    const Eval e = eval(x);
    if (distance(e.p1, e.p2) <= tol_.length) return x;
    return std::nullopt;
  }

  // Corrector: S1 = S2 and (S1 - P) . t = 0. Returns nullopt if not converged.
  std::optional<X4> correct(X4 x, const Vec3& P, const Vec3& t) const {
    for (int it = 0; it < 12; ++it) {
      const Eval e = eval(x);
      const Vec3 r = e.p1 - e.p2;
      const double g = dot(e.p1 - P, t);
      std::array<std::array<double, 4>, 4> A{};
      const std::array<Vec3, 4> J{e.su1, e.sv1, -e.su2, -e.sv2};
      for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t k = 0; k < 4; ++k) A[i][k] = comp(J[k], i);
      A[3] = {dot(t, e.su1), dot(t, e.sv1), 0.0, 0.0};
      std::array<double, 4> y{-r.x, -r.y, -r.z, -g};
      if (!solve_n<4>(A, y)) return std::nullopt;
      double step = 0.0;
      for (std::size_t k = 0; k < 4; ++k) {
        step += std::abs(y[k]) * norm(J[k]);
        x[k] += y[k];
      }
      x = clamp(x);
      if (step <= 1e-3 * tol_.length) break;
    }
    const Eval e = eval(x);
    if (distance(e.p1, e.p2) <= tol_.length) return x;
    return std::nullopt;
  }

  // Solve S1 = S2 with variable k fixed (on a domain boundary): 3 equations, 3 unknowns.
  std::optional<X4> solve_on_bound(X4 x, std::size_t k) const {
    for (int it = 0; it < 30; ++it) {
      const Eval e = eval(x);
      const Vec3 r = e.p1 - e.p2;
      const std::array<Vec3, 4> J{e.su1, e.sv1, -e.su2, -e.sv2};
      std::array<std::array<double, 3>, 3> A{};
      std::size_t col = 0;
      std::array<std::size_t, 3> idx{};
      for (std::size_t q = 0; q < 4; ++q) {
        if (q == k) continue;
        idx[col] = q;
        for (std::size_t i = 0; i < 3; ++i) A[i][col] = comp(J[q], i);
        ++col;
      }
      std::array<double, 3> y{-r.x, -r.y, -r.z};
      if (!solve_n<3>(A, y)) return std::nullopt;
      double step = 0.0;
      for (std::size_t c = 0; c < 3; ++c) {
        x[idx[c]] += y[c];
        step += std::abs(y[c]) * norm(J[idx[c]]);
      }
      x = clamp(x);
      if (step <= 1e-3 * tol_.length) break;
    }
    const Eval e = eval(x);
    if (distance(e.p1, e.p2) <= tol_.length) return x;
    return std::nullopt;
  }

  // ---------------------------------------------------------------- seeds
  void add_seed(const std::optional<X4>& x, std::vector<X4>& seeds) const {
    if (x) seeds.push_back(*x);
  }

  void collect_seeds(const NurbsSurface& p, const NurbsSurface& q, int depth, std::vector<X4>& seeds) {
    if (!overlap(box(p.control_points()), box(q.control_points()), tol_.length)) return;
    const Cone cp = normal_cone(p), cq = normal_cone(q);
    const double th = std::acos(std::clamp(dot(cp.axis, cq.axis), -1.0, 1.0));
    const double spread = cp.half + cq.half;
    const bool no_loops = cp.half < pi / 2 && cq.half < pi / 2 && th > spread && pi - th > spread;
    const double dp = distance(box(p.control_points()).lo, box(p.control_points()).hi);
    const double dq = distance(box(q.control_points()).lo, box(q.control_points()).hi);
    if (no_loops) {
      // every intersection curve in this pair crosses the boundary of p or q
      boundary_seeds(p, q, false, seeds);
      boundary_seeds(q, p, true, seeds);
      return;
    }
    if (depth >= 24 || std::max(dp, dq) <= 1e-4 * size_) {
      // tangential region: Newton seeds from a few samples
      const auto up = p.domain_u(), vp = p.domain_v(), uq = q.domain_u(), vq = q.domain_v();
      X4 best{};
      double bd = std::numeric_limits<double>::max();
      for (int i = 0; i <= 2; ++i)
        for (int j = 0; j <= 2; ++j) {
          const double u1 = up[0] + (up[1] - up[0]) * i / 2.0, v1 = vp[0] + (vp[1] - vp[0]) * j / 2.0;
          const Vec3 P = a_.evaluate(u1, v1);
          for (int k = 0; k <= 2; ++k)
            for (int l = 0; l <= 2; ++l) {
              const double u2 = uq[0] + (uq[1] - uq[0]) * k / 2.0, v2 = vq[0] + (vq[1] - vq[0]) * l / 2.0;
              const double d = distance(P, b_.evaluate(u2, v2));
              if (d < bd) { bd = d; best = {u1, v1, u2, v2}; }
            }
        }
      add_seed(refine(best), seeds);
      return;
    }
    if (dp >= dq) {
      const auto [l, r] = split(p);
      collect_seeds(l, q, depth + 1, seeds);
      collect_seeds(r, q, depth + 1, seeds);
    } else {
      const auto [l, r] = split(q);
      collect_seeds(p, l, depth + 1, seeds);
      collect_seeds(p, r, depth + 1, seeds);
    }
  }

  static std::pair<NurbsSurface, NurbsSurface> split(const NurbsSurface& s) {
    double ext_u = 0.0, ext_v = 0.0;
    for (std::size_t j = 0; j < s.num_v(); ++j)
      ext_u = std::max(ext_u, distance(s.control_point(0, j), s.control_point(s.num_u() - 1, j)));
    for (std::size_t i = 0; i < s.num_u(); ++i)
      ext_v = std::max(ext_v, distance(s.control_point(i, 0), s.control_point(i, s.num_v() - 1)));
    const auto ur = s.domain_u(), vr = s.domain_v();
    return ext_u >= ext_v ? s.split_u(0.5 * (ur[0] + ur[1])) : s.split_v(0.5 * (vr[0] + vr[1]));
  }

  // Boundary iso-curves of `edge_patch` intersected with `other`; `swapped` if edge_patch is a
  // piece of the second surface.
  void boundary_seeds(const NurbsSurface& edge_patch, const NurbsSurface& other, bool swapped,
                      std::vector<X4>& seeds) const {
    const auto ur = edge_patch.domain_u(), vr = edge_patch.domain_v();
    struct Edge {
      NurbsCurve c;
      bool along_u;  // curve parameter is u (v fixed)
      double fixed;
    };
    const std::array<Edge, 4> edges{Edge{edge_patch.iso_curve_v(vr[0]), true, vr[0]},
                                    Edge{edge_patch.iso_curve_v(vr[1]), true, vr[1]},
                                    Edge{edge_patch.iso_curve_u(ur[0]), false, ur[0]},
                                    Edge{edge_patch.iso_curve_u(ur[1]), false, ur[1]}};
    for (const auto& e : edges)
      for (const auto& h : intersect(e.c, other, tol_)) {
        const double eu = e.along_u ? h.t : e.fixed, ev = e.along_u ? e.fixed : h.t;
        const X4 x = swapped ? X4{h.u, h.v, eu, ev} : X4{eu, ev, h.u, h.v};
        add_seed(refine(x), seeds);
      }
  }

  // ---------------------------------------------------------------- tracing
  // Parameters of p and q agree, measured in model units (|du| |S_u| + |dv| |S_v| on both
  // surfaces). On closed surfaces the same 3D point occurs on both sides of a seam (large
  // parameter distance), while at a pole the u parameter is arbitrary (|S_u| = 0).
  bool params_near(const SurfaceSurfacePoint& p, const SurfaceSurfacePoint& q) const {
    const Eval e = eval({p.u1, p.v1, p.u2, p.v2});
    const double d1 = std::abs(p.u1 - q.u1) * norm(e.su1) + std::abs(p.v1 - q.v1) * norm(e.sv1);
    const double d2 = std::abs(p.u2 - q.u2) * norm(e.su2) + std::abs(p.v2 - q.v2) * norm(e.sv2);
    return std::max(d1, d2) <= 0.05 * size_;
  }

  bool near_traced(const SurfaceSurfacePoint& p) const {
    const double eps = 2.0 * chord_ + 10.0 * tol_.length;
    for (const auto& c : result_.curves)
      for (std::size_t i = 0; i + 1 < c.points.size(); ++i)
        if (segment_distance(p.point, c.points[i].point, c.points[i + 1].point) <= eps &&
            (params_near(p, c.points[i]) || params_near(p, c.points[i + 1])))
          return true;
    for (const auto& q : result_.points)
      if (distance(p.point, q.point) <= eps && params_near(p, q)) return true;
    return false;
  }

  SurfaceSurfacePoint make_point(const X4& x) const {
    const Eval e = eval(x);
    return {0.5 * (e.p1 + e.p2), x[0], x[1], x[2], x[3]};
  }

  // Marches from x0 in direction `sign` * tangent. Returns the points after x0 and whether the
  // curve closed back onto `start`.
  // A traced point; `brk` marks exact knot-line crossings (break points for the fit).
  struct Node {
    X4 x;
    bool brk;
  };

  // Points on the knot lines crossed between x and y (in order), appended to pts.
  void insert_knot_crossings(const X4& x, const X4& y, std::vector<Node>& pts) const {
    X4 from = x;
    for (int guard = 0; guard < 16; ++guard) {
      double fmin = 1.0;
      std::size_t kk = 4;
      double kappa = 0.0;
      for (std::size_t k = 0; k < 4; ++k)
        for (double kn : knots_[k]) {
          if ((from[k] - kn) * (y[k] - kn) >= 0.0) continue;  // not strictly crossed
          const double eps = 1e-9 * (hi_[k] - lo_[k]);
          if (std::abs(from[k] - kn) <= eps && std::abs(y[k] - kn) <= eps) continue;  // runs along it
          const double f = (kn - from[k]) / (y[k] - from[k]);
          if (f < fmin) { fmin = f; kk = k; kappa = kn; }
        }
      if (kk == 4) return;
      X4 guess;
      for (std::size_t k = 0; k < 4; ++k) guess[k] = from[k] + fmin * (y[k] - from[k]);
      guess[kk] = kappa;
      const auto xc = solve_on_bound(guess, kk);
      if (!xc) return;
      const Vec3 q = make_point(*xc).point;
      if (distance(q, make_point(from).point) > tol_.length && distance(q, make_point(y).point) > tol_.length)
        pts.push_back({*xc, true});
      from = *xc;
      from[kk] = kappa;
    }
  }

  std::vector<Node> march(const X4& x0, double sign, const Vec3& start, bool& closed) const {
    const X4 x_start = x0;
    std::vector<Node> pts;
    closed = false;
    X4 x = x0;
    auto t0 = tangent(x);
    if (!t0) return pts;
    Vec3 t = sign * *t0;
    const Vec3 t_start = t;
    double h = std::min(hmax_, 20.0 * chord_ + hmin_);
    double travelled = 0.0;
    for (int steps = 0; steps < 200000; ++steps) {
      const Eval e = eval(x);
      const Vec3 P = 0.5 * (e.p1 + e.p2);
      const auto d1 = param_dir(e.su1, e.sv1, t), d2 = param_dir(e.su2, e.sv2, t);
      const X4 dx{d1[0], d1[1], d2[0], d2[1]};
      bool accepted = false;
      while (h >= hmin_) {
        // predictor, clipped at the first domain boundary
        double f = 1.0;
        std::size_t hit = 4;
        for (std::size_t k = 0; k < 4; ++k) {
          const double xn = x[k] + h * dx[k];
          // leaving the domain by more than rounding noise (a curve may run along the boundary)
          const double eps = 1e-9 * (hi_[k] - lo_[k]);
          if (xn < lo_[k] - eps && dx[k] < 0.0) {
            const double fk = (lo_[k] - x[k]) / (h * dx[k]);
            if (fk < f) { f = fk; hit = k; }
          } else if (xn > hi_[k] + eps && dx[k] > 0.0) {
            const double fk = (hi_[k] - x[k]) / (h * dx[k]);
            if (fk < f) { f = fk; hit = k; }
          }
        }
        X4 xp;
        for (std::size_t k = 0; k < 4; ++k) xp[k] = x[k] + f * h * dx[k];
        if (hit < 4) {  // end on the boundary
          xp[hit] = dx[hit] < 0.0 ? lo_[hit] : hi_[hit];
          if (const auto xb = solve_on_bound(xp, hit)) {
            if (distance(make_point(*xb).point, P) > tol_.length) {
              insert_knot_crossings(x, *xb, pts);
              pts.push_back({*xb, false});
            }
            return pts;
          }
          h *= 0.5;
          continue;
        }
        const auto xc = correct(clamp(xp), P + h * t, t);
        if (!xc) { h *= 0.5; continue; }
        const Vec3 Pn = make_point(*xc).point;
        const double len = distance(Pn, P);
        auto tn = tangent(*xc);
        if (!tn || len < 0.3 * h || len > 2.0 * h) { h *= 0.5; continue; }
        if (dot(*tn, t) < 0.0) *tn = -*tn;
        const double theta = std::acos(std::clamp(dot(*tn, t), -1.0, 1.0));
        if (theta > 0.2 || len * theta / 8.0 > chord_) { h *= 0.5; continue; }
        // parameter-space control: near poles or strongly varying parametrisations the p-curves
        // bend much more than the 3D curve. Compare the corrected midpoint with the linear
        // parameter midpoint (in model units via |S_u|, |S_v|).
        {
          X4 xm;
          for (std::size_t k = 0; k < 4; ++k) xm[k] = 0.5 * (x[k] + (*xc)[k]);
          const Vec3 tm = normalized(t + *tn);
          if (const auto xs = correct(xm, 0.5 * (P + Pn), tm)) {
            const Eval em = eval(xm);
            const double dev1 = std::abs((*xs)[0] - xm[0]) * norm(em.su1) + std::abs((*xs)[1] - xm[1]) * norm(em.sv1);
            const double dev2 = std::abs((*xs)[2] - xm[2]) * norm(em.su2) + std::abs((*xs)[3] - xm[3]) * norm(em.sv2);
            if (std::max(dev1, dev2) > chord_) { h *= 0.5; continue; }
          }
        }
        // accept; first insert exact points where the step crosses knot lines of either surface
        // (the parametrisation is only C1/C2 there, so the p-curves are fitted piecewise)
        insert_knot_crossings(x, *xc, pts);
        x = *xc;
        t = *tn;
        travelled += len;
        pts.push_back({x, false});
        accepted = true;
        // closing the loop? (also in parameter space: on a closed surface the seam returns to
        // the same 3D point with other parameters)
        bool params_close = true;
        for (std::size_t k = 0; k < 4; ++k)
          if (std::abs(x[k] - x_start[k]) > 0.25 * (hi_[k] - lo_[k])) params_close = false;
        if (travelled > 3.0 * h && params_close && distance(Pn, start) < 1.2 * h && dot(t, t_start) > 0.0) {
          closed = true;
          return pts;
        }
        const double grow = theta > 0.0 ? std::sqrt(8.0 * chord_ / (len * theta)) : 2.0;
        h = std::clamp(h * std::clamp(grow, 0.5, 2.0), hmin_, hmax_);
        break;
      }
      if (!accepted) return pts;  // step size collapsed (singular point)
    }
    return pts;
  }

  void trace_from(const X4& seed) {
    const SurfaceSurfacePoint sp = make_point(seed);
    if (near_traced(sp)) return;
    if (!tangent(seed)) {  // tangential contact
      result_.points.push_back(sp);
      return;
    }
    bool closed = false, dummy = false;
    const auto fwd = march(seed, 1.0, sp.point, closed);
    std::vector<Node> bwd;
    if (!closed) bwd = march(seed, -1.0, sp.point, dummy);
    std::vector<Node> all(bwd.rbegin(), bwd.rend());
    all.push_back({seed, false});
    all.insert(all.end(), fwd.begin(), fwd.end());
    if (all.size() < 2) {
      result_.points.push_back(sp);
      return;
    }
    std::vector<SurfaceSurfacePoint> pts;
    std::vector<bool> brk;
    for (const auto& n : all) {
      pts.push_back(make_point(n.x));
      brk.push_back(n.brk);
    }
    if (closed) {
      pts.push_back(pts.front());
      brk.push_back(false);
    }
    result_.curves.push_back(fit(pts, brk, closed));
  }

  // Interpolates the 3D points and both parameter-space polylines with identical parameters,
  // piecewise between the knot-line crossings (where the surface parametrisations have reduced
  // continuity, a global cubic would converge only like h^2); the pieces are joined (C0).
  // Short pieces (closely spaced knot lines) would be interpolated with degree < 3; they get
  // additional exact points on the intersection curve (corrected midpoints of the longest gap).
  // Exact intersection point between two traced points (corrected parameter midpoint).
  std::optional<SurfaceSurfacePoint> midpoint(const SurfaceSurfacePoint& p, const SurfaceSurfacePoint& q) const {
    const X4 xm{0.5 * (p.u1 + q.u1), 0.5 * (p.v1 + q.v1), 0.5 * (p.u2 + q.u2), 0.5 * (p.v2 + q.v2)};
    const Vec3 chord = q.point - p.point;
    if (norm(chord) <= 0.0) return std::nullopt;
    const auto xs = correct(xm, 0.5 * (p.point + q.point), chord / norm(chord));
    if (!xs) return std::nullopt;
    return make_point(*xs);
  }

  struct PieceFit {
    NurbsCurve c, c1, c2;
  };
  static PieceFit fit_piece(const std::vector<SurfaceSurfacePoint>& piece, double t0) {
    std::vector<Vec3> P, Q1, Q2;
    for (const auto& q : piece) {
      P.push_back(q.point);
      Q1.push_back({q.u1, q.v1, 0.0});
      Q2.push_back({q.u2, q.v2, 0.0});
    }
    std::vector<double> t(P.size(), t0);  // chord length, continuing from the previous piece
    for (std::size_t k = 1; k < P.size(); ++k) t[k] = t[k - 1] + distance(P[k], P[k - 1]);
    const std::size_t deg = std::min<std::size_t>(3, P.size() - 1);
    return {interpolate_curve(P, deg, t), interpolate_curve(Q1, deg, t), interpolate_curve(Q2, deg, t)};
  }

  // Fits one piece; short pieces get extra points (degree 3), and exact points are inserted
  // between fit points until |S_i(pcurve_i(t)) - curve(t)| <= fit tolerance at all midpoints.
  PieceFit fit_verified(std::vector<SurfaceSurfacePoint>& piece, double t0) const {
    for (int iter = 0; iter < 40; ++iter) {
      if (piece.size() < 4) {  // longest gap first
        std::size_t g = 0;
        for (std::size_t k = 1; k + 1 < piece.size(); ++k)
          if (distance(piece[k + 1].point, piece[k].point) > distance(piece[g + 1].point, piece[g].point)) g = k;
        if (const auto m = midpoint(piece[g], piece[g + 1])) {
          piece.insert(piece.begin() + static_cast<std::ptrdiff_t>(g + 1), *m);
          continue;
        }
      }
      PieceFit f = fit_piece(piece, t0);
      std::vector<std::size_t> bad;
      double tk = t0;
      for (std::size_t k = 0; k + 1 < piece.size(); ++k) {
        const double tn = tk + distance(piece[k + 1].point, piece[k].point);
        const double tm = 0.5 * (tk + tn);
        const Vec3 q = f.c.evaluate(tm), w1 = f.c1.evaluate(tm), w2 = f.c2.evaluate(tm);
        const double e = std::max(distance(a_.evaluate(w1.x, w1.y), q), distance(b_.evaluate(w2.x, w2.y), q));
        if (e > fit_) bad.push_back(k);
        tk = tn;
      }
      if (bad.empty() || piece.size() > 20000) return f;
      bool inserted = false;
      for (auto it = bad.rbegin(); it != bad.rend(); ++it)
        if (const auto m = midpoint(piece[*it], piece[*it + 1])) {
          piece.insert(piece.begin() + static_cast<std::ptrdiff_t>(*it + 1), *m);
          inserted = true;
        }
      if (!inserted) return f;
    }
    return fit_piece(piece, t0);
  }

  IntersectionCurve fit(const std::vector<SurfaceSurfacePoint>& pts, const std::vector<bool>& brk,
                        bool closed) const {
    std::optional<NurbsCurve> c, c1, c2;
    std::vector<SurfaceSurfacePoint> all;
    std::size_t begin = 0;
    for (std::size_t i = 1; i < pts.size(); ++i) {
      if (!brk[i] && i + 1 != pts.size()) continue;
      std::vector<SurfaceSurfacePoint> piece(pts.begin() + static_cast<std::ptrdiff_t>(begin),
                                             pts.begin() + static_cast<std::ptrdiff_t>(i + 1));
      PieceFit f = fit_verified(piece, c ? c->domain_end() : 0.0);
      c = c ? join(*c, f.c, 1e-9) : f.c;
      c1 = c1 ? join(*c1, f.c1, 1e-9) : f.c1;
      c2 = c2 ? join(*c2, f.c2, 1e-9) : f.c2;
      all.insert(all.end(), piece.begin() + (all.empty() ? 0 : 1), piece.end());
      begin = i;
    }
    return {std::move(all), closed, std::move(*c), std::move(*c1), std::move(*c2)};
  }

  const NurbsSurface& a_;
  const NurbsSurface& b_;
  Tolerance tol_;
  double size_{1.0}, chord_{1e-4}, fit_{1e-7}, hmax_{1.0}, hmin_{1e-9};
  X4 lo_{}, hi_{};
  std::array<std::vector<double>, 4> knots_;  // interior knot values per parameter
  SurfaceIntersection result_;
};

}  // namespace

SurfaceIntersection intersect(const NurbsSurface& a, const NurbsSurface& b,
                              const SurfaceIntersectionOptions& options) {
  return Ssi(a, b, options).run();
}

}  // namespace nf
