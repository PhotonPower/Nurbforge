#include "nurbforge/analysis/intersect.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <vector>
#include <numbers>

#include "nurbforge/analysis/closest_point.hpp"
#include "nurbforge/core/error.hpp"

namespace nf {

namespace {

constexpr double kNoise = 1.0 + 1e-12;  // backtracking tolerance for rounding noise
constexpr int kMaxDepth = 60;
constexpr double pi = std::numbers::pi;

struct Box3 {
  Vec3 lo, hi;
};

Box3 box(const std::vector<Vec3>& pts, double inflate) {
  Box3 b{pts.front(), pts.front()};
  for (const auto& p : pts) {
    b.lo = {std::min(b.lo.x, p.x), std::min(b.lo.y, p.y), std::min(b.lo.z, p.z)};
    b.hi = {std::max(b.hi.x, p.x), std::max(b.hi.y, p.y), std::max(b.hi.z, p.z)};
  }
  b.lo = b.lo - Vec3{inflate, inflate, inflate};
  b.hi = b.hi + Vec3{inflate, inflate, inflate};
  return b;
}

bool overlap(const Box3& a, const Box3& b) {
  return a.lo.x <= b.hi.x && b.lo.x <= a.hi.x && a.lo.y <= b.hi.y && b.lo.y <= a.hi.y &&
         a.lo.z <= b.hi.z && b.lo.z <= a.hi.z;
}

double diagonal(const std::vector<Vec3>& pts) {
  const Box3 b = box(pts, 0.0);
  return distance(b.lo, b.hi);
}

// Cone containing all directions of the given edge vectors: axis and half angle (pi if none).
struct Cone {
  Vec3 axis;
  double half{pi};
};

Cone cone(const std::vector<Vec3>& edges, const Vec3& axis_guess) {
  const double la = norm(axis_guess);
  if (la <= 0.0) return {};
  Cone c{axis_guess / la, 0.0};
  for (const auto& e : edges) {
    const double le = norm(e);
    if (le == 0.0) continue;
    c.half = std::max(c.half, std::acos(std::clamp(dot(e, c.axis) / le, -1.0, 1.0)));
  }
  return c;
}

double angle(const Vec3& a, const Vec3& b) { return std::acos(std::clamp(dot(a, b), -1.0, 1.0)); }

Cone curve_cone(const NurbsCurve& c) {
  const auto& P = c.control_points();
  std::vector<Vec3> e;
  for (std::size_t i = 1; i < P.size(); ++i) e.push_back(P[i] - P[i - 1]);
  return cone(e, P.back() - P.front());
}

// Cone of the normals of a patch (conservative estimate from the u- and v-edge cones).
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
  const Cone cu = cone(eu, du), cv = cone(ev, dv);
  const Vec3 n = cross(du, dv);
  if (norm(n) == 0.0 || cu.half >= pi / 2 || cv.half >= pi / 2) return {};
  const double sin_psi = norm(cross(cu.axis, cv.axis));
  return {n / norm(n), (cu.half + cv.half) / std::max(sin_psi, 0.1)};
}

// Small dense solve (N <= 3) with partial pivoting; returns false if singular.
template <std::size_t N>
bool solve(std::array<std::array<double, N>, N> A, std::array<double, N>& b, std::size_t n) {
  for (std::size_t k = 0; k < n; ++k) {
    std::size_t p = k;
    for (std::size_t i = k + 1; i < n; ++i)
      if (std::abs(A[i][k]) > std::abs(A[p][k])) p = i;
    if (std::abs(A[p][k]) < 1e-300) return false;
    std::swap(A[k], A[p]);
    std::swap(b[k], b[p]);
    for (std::size_t i = k + 1; i < n; ++i) {
      const double f = A[i][k] / A[k][k];
      for (std::size_t j = k; j < n; ++j) A[i][j] -= f * A[k][j];
      b[i] -= f * b[k];
    }
  }
  for (std::size_t i = n; i-- > 0;) {
    for (std::size_t j = i + 1; j < n; ++j) b[i] -= A[i][j] * b[j];
    b[i] /= A[i][i];
  }
  return true;
}

// Box-constrained Gauss-Newton for r(x) = 0 with r: R^N -> R^3; `eval` returns r and the
// Jacobian columns. Projected (active set) steps, backtracking on |r|^2.
template <std::size_t N, class Eval>
std::array<double, N> gauss_newton(std::array<double, N> x, const std::array<double, N>& lo,
                                   const std::array<double, N>& hi, const Eval& eval) {
  for (int it = 0; it < 80; ++it) {
    std::array<Vec3, N> J;
    const Vec3 r = eval(x, &J);
    const double r2 = dot(r, r);
    if (r2 == 0.0) break;
    std::array<double, N> g{};
    for (std::size_t i = 0; i < N; ++i) g[i] = dot(J[i], r);
    // active set: variables on a bound whose descent direction points outwards
    std::array<std::size_t, N> free{};
    std::size_t nf = 0;
    for (std::size_t i = 0; i < N; ++i)
      if (!((x[i] <= lo[i] && g[i] > 0.0) || (x[i] >= hi[i] && g[i] < 0.0))) free[nf++] = i;
    if (nf == 0) break;
    std::array<std::array<double, N>, N> A{};
    std::array<double, N> rhs{};
    double trace = 0.0;
    for (std::size_t a = 0; a < nf; ++a) {
      for (std::size_t b = 0; b < nf; ++b) A[a][b] = dot(J[free[a]], J[free[b]]);
      rhs[a] = -g[free[a]];
      trace += A[a][a];
    }
    for (std::size_t a = 0; a < nf; ++a) A[a][a] += 1e-14 * trace;  // tangential / singular cases
    if (!solve<N>(A, rhs, nf)) break;
    std::array<double, N> dx{};
    for (std::size_t a = 0; a < nf; ++a) dx[free[a]] = rhs[a];
    double t = 1.0;
    std::array<double, N> xn = x;
    while (true) {
      for (std::size_t i = 0; i < N; ++i) xn[i] = std::clamp(x[i] + t * dx[i], lo[i], hi[i]);
      const Vec3 rn = eval(xn, nullptr);
      if (dot(rn, rn) <= r2 * kNoise || t < 1e-10) break;
      t *= 0.5;
    }
    double step = 0.0;
    for (std::size_t i = 0; i < N; ++i) step += std::abs(xn[i] - x[i]) * norm(J[i]);
    x = xn;
    if (step <= 1e-15 * (1.0 + std::sqrt(r2)) + 1e-300) break;
  }
  return x;
}

// Number of grid start points tried by the leaf Gauss-Newton solves.
constexpr std::size_t kStarts = 4;

// ------------------------------------------------------------------ curve / curve

struct CurveCurveSearch {
  const NurbsCurve& a;
  const NurbsCurve& b;
  const Tolerance& tol;
  double min_size;
  std::vector<CurveCurveIntersection> out;

  void leaf(const NurbsCurve& pa, const NurbsCurve& pb) {
    const double a0 = pa.domain_start(), a1 = pa.domain_end(), b0 = pb.domain_start(), b1 = pb.domain_end();
    std::vector<std::pair<double, std::array<double, 2>>> starts;
    for (int i = 0; i <= 4; ++i)
      for (int j = 0; j <= 4; ++j) {
        const double s = a0 + (a1 - a0) * i / 4.0, t = b0 + (b1 - b0) * j / 4.0;
        starts.push_back({distance(a.evaluate(s), b.evaluate(t)), {s, t}});
      }
    std::sort(starts.begin(), starts.end(), [](const auto& l, const auto& r) { return l.first < r.first; });
    const auto eval = [&](const std::array<double, 2>& p, std::array<Vec3, 2>* J) {
      if (!J) return a.evaluate(p[0]) - b.evaluate(p[1]);
      const auto da = a.derivatives(p[0], 1), db = b.derivatives(p[1], 1);
      (*J)[0] = da[1];
      (*J)[1] = -db[1];
      return da[0] - db[0];
    };
    // Gauss-Newton from the best grid points: from a single start the box-constrained iteration
    // can stall on a bound when the pieces are large (e.g. separated cones at the first level)
    for (std::size_t k = 0; k < std::min<std::size_t>(kStarts, starts.size()); ++k) {
      const auto x = gauss_newton<2>(starts[k].second, {a0, b0}, {a1, b1}, eval);
      const Vec3 pa_ = a.evaluate(x[0]), pb_ = b.evaluate(x[1]);
      const double d = distance(pa_, pb_);
      if (d <= tol.length) {
        out.push_back({x[0], x[1], 0.5 * (pa_ + pb_), d});
        return;
      }
    }
  }

  // Overlap: if the piece `p` lies on the whole curve `q` (5 samples within tolerance, local
  // projections started at `guess` and continued from sample to sample), its end points are
  // reported and the pair is not subdivided further (otherwise coincident curves would be
  // subdivided down to the size limit everywhere). A failed local projection only means more
  // subdivision, never a wrong result.
  bool covered(const NurbsCurve& p, const NurbsCurve& q, bool p_is_b, double guess) {
    std::array<CurvePoint, 5> cp{};
    for (int i = 0; i <= 4; ++i) {
      const double t = p.domain_start() + (p.domain_end() - p.domain_start()) * i / 4.0;
      cp[static_cast<std::size_t>(i)] = closest_point_local(q, p.evaluate(t), guess, tol);
      guess = cp[static_cast<std::size_t>(i)].u;
      if (cp[static_cast<std::size_t>(i)].distance > tol.length) return false;
    }
    for (int i : {0, 4}) {
      const double tp = i == 0 ? p.domain_start() : p.domain_end();
      const auto& c = cp[static_cast<std::size_t>(i)];
      const double s_ = p_is_b ? c.u : tp, t_ = p_is_b ? tp : c.u;
      out.push_back({s_, t_, 0.5 * (a.evaluate(s_) + b.evaluate(t_)), c.distance});
    }
    return true;
  }

  void process(const NurbsCurve& pa, const NurbsCurve& pb, int depth) {
    const auto& A = pa.control_points();
    const auto& B = pb.control_points();
    if (!overlap(box(A, tol.length), box(B, tol.length))) return;
    const double da = diagonal(A), db = diagonal(B);
    const Cone ca = curve_cone(pa), cb = curve_cone(pb);
    const double th = angle(ca.axis, cb.axis), spread = ca.half + cb.half;
    const bool separated = ca.half < pi / 2 && cb.half < pi / 2 && th > spread && pi - th > spread;
    if (separated || depth >= kMaxDepth || (da <= min_size && db <= min_size)) {
      leaf(pa, pb);
      return;
    }
    const double mid_a = 0.5 * (pa.domain_start() + pa.domain_end());
    const double mid_b = 0.5 * (pb.domain_start() + pb.domain_end());
    if (covered(pa, b, false, mid_b) || covered(pb, a, true, mid_a)) return;
    if (da >= db) {
      const auto [l, r] = pa.split(0.5 * (pa.domain_start() + pa.domain_end()));
      process(l, pb, depth + 1);
      process(r, pb, depth + 1);
    } else {
      const auto [l, r] = pb.split(0.5 * (pb.domain_start() + pb.domain_end()));
      process(pa, l, depth + 1);
      process(pa, r, depth + 1);
    }
  }
};

// ------------------------------------------------------------------ curve / surface

struct CurveSurfaceSearch {
  const NurbsCurve& c;
  const NurbsSurface& s;
  const Tolerance& tol;
  double min_size;
  std::vector<CurveSurfaceIntersection> out;

  void leaf(const NurbsCurve& pc, const NurbsSurface& ps) {
    const double t0 = pc.domain_start(), t1 = pc.domain_end();
    const auto ur = ps.domain_u(), vr = ps.domain_v();
    std::vector<std::pair<double, std::array<double, 3>>> starts;
    for (int i = 0; i <= 4; ++i) {
      const double t = t0 + (t1 - t0) * i / 4.0;
      const Vec3 q = c.evaluate(t);
      for (int j = 0; j <= 4; ++j)
        for (int k = 0; k <= 4; ++k) {
          const double u = ur[0] + (ur[1] - ur[0]) * j / 4.0, v = vr[0] + (vr[1] - vr[0]) * k / 4.0;
          starts.push_back({distance(q, s.evaluate(u, v)), {t, u, v}});
        }
    }
    std::sort(starts.begin(), starts.end(), [](const auto& l, const auto& r) { return l.first < r.first; });
    const auto eval = [&](const std::array<double, 3>& p, std::array<Vec3, 3>* J) {
      if (!J) return c.evaluate(p[0]) - s.evaluate(p[1], p[2]);
      const auto dc = c.derivatives(p[0], 1);
      const auto ds = s.first_derivatives(p[1], p[2]);
      (*J)[0] = dc[1];
      (*J)[1] = -ds[1];
      (*J)[2] = -ds[2];
      return dc[0] - ds[0];
    };
    for (std::size_t k = 0; k < std::min<std::size_t>(kStarts, starts.size()); ++k) {
      const auto x = gauss_newton<3>(starts[k].second, {t0, ur[0], vr[0]}, {t1, ur[1], vr[1]}, eval);
      const Vec3 pc_ = c.evaluate(x[0]), ps_ = s.evaluate(x[1], x[2]);
      const double d = distance(pc_, ps_);
      if (d <= tol.length) {
        out.push_back({x[0], x[1], x[2], 0.5 * (pc_ + ps_), d});
        return;
      }
    }
  }

  // Curve piece lying in the surface (5 samples within tolerance, local projections from the
  // piece's centre): report its end points instead of subdividing a coincident configuration
  // down to the size limit.
  bool covered(const NurbsCurve& pc, const NurbsSurface& ps) {
    const auto ur = ps.domain_u(), vr = ps.domain_v();
    double gu = 0.5 * (ur[0] + ur[1]), gv = 0.5 * (vr[0] + vr[1]);
    std::array<SurfacePoint, 5> sp{};
    const double band = 10.0 * min_size;
    bool within_tol = true;
    for (int i = 0; i <= 4; ++i) {
      const double t = pc.domain_start() + (pc.domain_end() - pc.domain_start()) * i / 4.0;
      sp[static_cast<std::size_t>(i)] = closest_point_local(s, pc.evaluate(t), gu, gv, tol);
      gu = sp[static_cast<std::size_t>(i)].u;
      gv = sp[static_cast<std::size_t>(i)].v;
      if (sp[static_cast<std::size_t>(i)].distance > band) return false;
      if (sp[static_cast<std::size_t>(i)].distance > tol.length) within_tol = false;
    }
    if (!within_tol) return near_parallel(pc, sp[0].u, sp[0].v);
    for (int i : {0, 4}) {
      const double t = i == 0 ? pc.domain_start() : pc.domain_end();
      const auto& q = sp[static_cast<std::size_t>(i)];
      out.push_back({t, q.u, q.v, 0.5 * (c.evaluate(t) + s.evaluate(q.u, q.v)), q.distance});
    }
    return true;
  }

  // Curve piece running close to the surface without lying in it (e.g. a line just off a
  // tilted plane, or the boundary of a patch next to a line of tangential contact): the boxes
  // overlap all along it, and subdivision would go down to the size limit everywhere. Instead
  // the signed distance is scanned along the piece; sign changes are refined by bisection,
  // minima within the tolerance are contacts. Returns false (subdivide) if a projection fails.
  bool near_parallel(const NurbsCurve& pc, double u0, double v0) {
    constexpr int n = 64;
    const double t0 = pc.domain_start(), t1 = pc.domain_end();
    struct Sample {
      double t, d;
      SurfacePoint f;
    };
    double gu = u0, gv = v0;
    const auto signed_distance = [&](double t, Sample& out_s) {
      const Vec3 P = pc.evaluate(t);
      const SurfacePoint f = closest_point_local(s, P, gu, gv, tol);
      if (f.distance > 20.0 * min_size) return false;
      Vec3 nrm;
      try {
        nrm = s.normal(f.u, f.v);
      } catch (const GeometryError&) {
        return false;
      }
      gu = f.u;
      gv = f.v;
      out_s = {t, dot(P - f.point, nrm), f};
      return true;
    };
    std::vector<Sample> smp(n + 1);
    for (int i = 0; i <= n; ++i)
      if (!signed_distance(t0 + (t1 - t0) * i / n, smp[static_cast<std::size_t>(i)])) return false;
    const auto push = [&](const Sample& x) {
      const Vec3 P = c.evaluate(x.t);
      out.push_back({x.t, x.f.u, x.f.v, 0.5 * (P + x.f.point), std::abs(x.d)});
    };
    for (std::size_t i = 0; i < smp.size(); ++i) {
      const Sample& x = smp[i];
      if (std::abs(x.d) <= tol.length) {
        push(x);
        continue;
      }
      if (i + 1 < smp.size() && x.d * smp[i + 1].d < 0.0 && std::abs(smp[i + 1].d) > tol.length) {
        Sample lo = x, hi = smp[i + 1], mid{};
        for (int it = 0; it < 60; ++it) {
          gu = lo.f.u;
          gv = lo.f.v;
          if (!signed_distance(0.5 * (lo.t + hi.t), mid)) return false;
          if (mid.d == 0.0 || hi.t - lo.t <= 1e-15 * (t1 - t0)) break;
          (mid.d * lo.d > 0.0 ? lo : hi) = mid;
        }
        push(mid);
      }
    }
    return true;
  }

  void process(const NurbsCurve& pc, const NurbsSurface& ps, int depth) {
    const auto& C = pc.control_points();
    const auto& S = ps.control_points();
    if (!overlap(box(C, tol.length), box(S, tol.length))) return;
    const double dc = diagonal(C), dsz = diagonal(S);
    const Cone cc = curve_cone(pc), cn = normal_cone(ps);
    bool transversal = false;
    if (cc.half < pi / 2 && cn.half < pi / 2) {
      const double phi = angle(cc.axis, cn.axis);
      transversal = std::min(phi, pi - phi) + cc.half < pi / 2 - cn.half;
    }
    if (transversal || depth >= kMaxDepth || (dc <= min_size && dsz <= min_size)) {
      leaf(pc, ps);
      return;
    }
    if (covered(pc, ps)) return;
    if (dc >= dsz) {
      const auto [l, r] = pc.split(0.5 * (pc.domain_start() + pc.domain_end()));
      process(l, ps, depth + 1);
      process(r, ps, depth + 1);
      return;
    }
    // split the patch along its larger extent
    double ext_u = 0.0, ext_v = 0.0;
    for (std::size_t j = 0; j < ps.num_v(); ++j)
      ext_u = std::max(ext_u, distance(ps.control_point(0, j), ps.control_point(ps.num_u() - 1, j)));
    for (std::size_t i = 0; i < ps.num_u(); ++i)
      ext_v = std::max(ext_v, distance(ps.control_point(i, 0), ps.control_point(i, ps.num_v() - 1)));
    const auto ur = ps.domain_u(), vr = ps.domain_v();
    const auto [l, r] = ext_u >= ext_v ? ps.split_u(0.5 * (ur[0] + ur[1])) : ps.split_v(0.5 * (vr[0] + vr[1]));
    process(pc, l, depth + 1);
    process(pc, r, depth + 1);
  }
};

// Merges candidates of the same contact: nearby points, or consecutive (in the curve parameter)
// candidates whose parameter midpoint is still within tolerance – a tangential contact keeps
// the distance below tol over a short stretch (d ~ kappa dt^2 / 2), where Gauss-Newton may stop
// anywhere. `mid_gap(a, b)` returns the gap at the midpoint between two candidates.
template <class T, class Key, class MidGap>
std::vector<T> deduplicate(std::vector<T> v, double merge_dist, double tol, Key key, MidGap mid_gap) {
  std::sort(v.begin(), v.end(), [&](const T& x, const T& y) { return key(x) < key(y); });
  std::vector<T> out;
  for (const auto& p : v) {
    bool merged = false;
    for (auto& q : out)
      if (distance(p.point, q.point) <= merge_dist) {
        if (p.distance < q.distance) q = p;
        merged = true;
        break;
      }
    if (!merged && !out.empty() && mid_gap(out.back(), p) <= tol) {
      if (p.distance < out.back().distance) out.back() = p;
      merged = true;
    }
    if (!merged) out.push_back(p);
  }
  return out;
}

}  // namespace

std::vector<CurveCurveIntersection> intersect(const NurbsCurve& a, const NurbsCurve& b, const Tolerance& tol) {
  std::vector<Vec3> all = a.control_points();
  all.insert(all.end(), b.control_points().begin(), b.control_points().end());
  const double scale = diagonal(all);
  CurveCurveSearch search{a, b, tol, std::max(tol.length, 1e-6 * scale), {}};
  for (const auto& pa : a.decompose_bezier())
    for (const auto& pb : b.decompose_bezier()) search.process(pa, pb, 0);
  const auto mid_gap = [&](const CurveCurveIntersection& x, const CurveCurveIntersection& y) {
    const Vec3 q = a.evaluate(0.5 * (x.s + y.s));
    return closest_point_local(b, q, 0.5 * (x.t + y.t), tol).distance;
  };
  return deduplicate(std::move(search.out), std::max(10 * tol.length, 1e-7 * scale), tol.length,
                     [](const CurveCurveIntersection& x) { return x.s; }, mid_gap);
}

std::vector<CurveSurfaceIntersection> intersect(const NurbsCurve& c, const NurbsSurface& s, const Tolerance& tol) {
  std::vector<Vec3> all = c.control_points();
  all.insert(all.end(), s.control_points().begin(), s.control_points().end());
  const double scale = diagonal(all);
  CurveSurfaceSearch search{c, s, tol, std::max(tol.length, 1e-6 * scale), {}};
  const auto patches = s.decompose_bezier();
  for (const auto& pc : c.decompose_bezier())
    for (const auto& ps : patches) search.process(pc, ps, 0);
  const auto mid_gap = [&](const CurveSurfaceIntersection& x, const CurveSurfaceIntersection& y) {
    const Vec3 q = c.evaluate(0.5 * (x.t + y.t));
    return closest_point_local(s, q, 0.5 * (x.u + y.u), 0.5 * (x.v + y.v), tol).distance;
  };
  return deduplicate(std::move(search.out), std::max(10 * tol.length, 1e-7 * scale), tol.length,
                     [](const CurveSurfaceIntersection& x) { return x.t; }, mid_gap);
}

}  // namespace nf
