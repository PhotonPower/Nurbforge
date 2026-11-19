#include "nurbforge/topology/edge_blend.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <iterator>
#include <limits>
#include <map>
#include <memory>
#include <numbers>
#include <optional>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "../modeling/geodesic.hpp"
#include "nurbforge/analysis/closest_point.hpp"
#include "nurbforge/analysis/intersect.hpp"
#include "nurbforge/analysis/properties.hpp"
#include "nurbforge/core/error.hpp"
#include "nurbforge/fitting/fit.hpp"
#include "nurbforge/geometry/primitives.hpp"
#include "nurbforge/modeling/construct.hpp"
#include "nurbforge/modeling/extend.hpp"
#include "nurbforge/modeling/fillet.hpp"
#include "nurbforge/modeling/project.hpp"
#include "nurbforge/nurbs/join.hpp"
#include "nurbforge/nurbs/knot_vector.hpp"

namespace nf {

namespace {

std::string num(double x) { return std::to_string(x); }

const TrimLoop& loop_of(const TrimmedSurface& f, std::size_t l) {
  return l == 0 ? f.outer() : f.holes()[l - 1];
}

Vec3 at_start(const NurbsCurve& c) { return c.evaluate(c.domain_start()); }
Vec3 at_end(const NurbsCurve& c) { return c.evaluate(c.domain_end()); }
Vec3 on(const NurbsSurface& s, const Vec3& uv) { return s.evaluate(uv.x, uv.y); }

// The curve with its end points moved exactly onto p0 / p1 (clamped curves interpolate their
// end control points); used to close loops whose pieces were computed separately.
NurbsCurve snapped(const NurbsCurve& c, const Vec3& p0, const Vec3& p1) {
  const NurbsCurve k = c.clamped();
  auto P = k.control_points();
  P.front() = p0;
  P.back() = p1;
  return {k.degree(), k.knots(), std::move(P), k.weights()};
}

// Opposite boundaries of the surface coincide (closed or periodic direction).
bool closed_in(const NurbsSurface& s, bool in_u) {
  const auto du = s.domain_u(), dv = s.domain_v();
  double size = 0.0;
  for (const auto& p : s.control_points())
    size = std::max(size, distance(p, s.control_points()[0]));
  for (int i = 0; i <= 4; ++i) {
    const double t = i / 4.0;
    const Vec3 a = in_u ? s.evaluate(du[0], dv[0] + t * (dv[1] - dv[0]))
                        : s.evaluate(du[0] + t * (du[1] - du[0]), dv[0]);
    const Vec3 b = in_u ? s.evaluate(du[1], dv[0] + t * (dv[1] - dv[0]))
                        : s.evaluate(du[0] + t * (du[1] - du[0]), dv[1]);
    if (distance(a, b) > 1e-9 * std::max(size, 1.0)) return false;
  }
  return true;
}

// The surface extended by `m` across every side that is not a seam (parameters kept): linearly
// (along the tangent planes) where the iso-curves across the side are polynomial - a natural
// extension of a fitted surface with short spans bends away quickly - and naturally otherwise
// (rational: conics continue exactly).
NurbsSurface extended(const NurbsSurface& s, double m) {
  NurbsSurface r = s;
  const auto grow = [&](SurfaceSide side) {
    try {
      r = extend(r, side, m, Extension::Linear);
      return;
    } catch (const GeometryError&) {
    }
    try {
      r = extend(r, side, m);
    } catch (const GeometryError&) {  // degenerate side (pole): leave it
    }
  };
  if (!closed_in(s, true)) {
    grow(SurfaceSide::UMin);
    grow(SurfaceSide::UMax);
  }
  if (!closed_in(s, false)) {
    grow(SurfaceSide::VMin);
    grow(SurfaceSide::VMax);
  }
  return r;
}

// The curve extended at one end like the surfaces in extended().
NurbsCurve extended_curve(const NurbsCurve& c, CurveEnd end, double m) {
  if (!c.is_rational()) return extend(c, end, m, Extension::Linear);
  return extend(c, end, m);
}

// Adaptive interpolation of exactly computed points: params T are refined (midpoints added)
// until every interpolated curve is within tol of the exact values at the midpoints. `eval`
// returns, for a parameter, the exact points of all curves; `check` measures the deviation of
// the interpolations from them at a parameter (in model units).
std::vector<NurbsCurve> fit_family(
    std::vector<double> T, std::size_t ncurves,
    const std::function<std::vector<Vec3>(double)>& eval,
    const std::function<double(const std::vector<NurbsCurve>&, double, const std::vector<Vec3>&)>&
        check,
    double tol, const std::string& who) {
  std::map<double, std::vector<Vec3>> cache;
  const auto exact = [&](double t) -> const std::vector<Vec3>& {
    auto it = cache.find(t);
    if (it == cache.end()) it = cache.emplace(t, eval(t)).first;
    return it->second;
  };
  for (int round = 0; round < 14; ++round) {
    std::vector<std::vector<Vec3>> pts(ncurves);
    for (const double t : T) {
      const auto& e = exact(t);
      for (std::size_t k = 0; k < ncurves; ++k) pts[k].push_back(e[k]);
    }
    const std::size_t deg = std::min<std::size_t>(3, T.size() - 1);
    std::vector<NurbsCurve> curves;
    for (std::size_t k = 0; k < ncurves; ++k) curves.push_back(interpolate_curve(pts[k], deg, T));
    std::vector<double> add;
    for (std::size_t i = 0; i + 1 < T.size(); ++i) {
      const double tm = 0.5 * (T[i] + T[i + 1]);
      if (check(curves, tm, exact(tm)) > tol) add.push_back(tm);
    }
    if (add.empty()) return curves;
    if (T.size() + add.size() > 4097) break;
    T.insert(T.end(), add.begin(), add.end());
    std::sort(T.begin(), T.end());
  }
  throw GeometryError(who + ": tolerance not reached (surfaces not smooth enough?)");
}

std::vector<double> uniform(double a, double b, int n) {
  std::vector<double> T;
  for (int i = 0; i <= n; ++i) T.push_back(i == n ? b : a + (b - a) * i / n);
  return T;
}

// ------------------------------------------------------------------ the blend between faces

// Blend surface B(u, v): u in [0, 1] across (u = 0 on face 1, u = 1 on face 2), v along the
// edge; k1 / k2 the 3D contact curves B(0, v) / B(1, v) and pc1 / pc2 their p-curves on the
// two face surfaces, all with the parameter v.
struct Blend {
  NurbsSurface surface;
  NurbsCurve k1, k2, pc1, pc2;
};

// A group of selected edges between the same two faces, resolved on the current shell.
struct Group {
  std::size_t f1{0}, f2{0};
  std::vector<std::size_t> edges;
  std::vector<std::size_t> ends;  // vertex indices at the ends of an open chain
};

struct Coedges {
  std::size_t loop{0};
  std::vector<std::size_t> idx;  // positions in the loop (unordered)
};

Coedges coedges_in(const Shell& sh, const Group& g, std::size_t face) {
  Coedges c;
  bool first = true;
  for (const std::size_t e : g.edges)
    for (const auto& u : sh.edges()[e].uses)
      if (u.face == face) {
        if (!first && u.loop != c.loop)
          throw GeometryError("edge blend: the edges lie on different boundary loops of a face");
        c.loop = u.loop;
        c.idx.push_back(u.index);
        first = false;
      }
  return c;
}

// Material side: a convex edge (dihedral angle inside the material below 180 degrees).
bool is_convex(const Shell& sh, const ShellEdge& e) {
  const auto& F1 = sh.faces()[e.uses[0].face];
  const auto& F2 = sh.faces()[e.uses[1].face];
  const NurbsCurve& pc = loop_of(F1, e.uses[0].loop).curves[e.uses[0].index];
  const double tm = 0.5 * (pc.domain_start() + pc.domain_end());
  const auto d = pc.derivatives(tm, 1);
  const auto D = F1.surface().first_derivatives(d[0].x, d[0].y);
  const Vec3 T = d[1].x * D[1] + d[1].y * D[2];
  const Vec3 into1 = cross(cross(D[1], D[2]), T);  // left of the loop: face interior
  const SurfacePoint q = closest_point(F2.surface(), D[0]);
  const Vec3 n2 = F2.surface().normal(q.u, q.v);
  return dot(into1, n2) < 0.0;
}

// Geodesic chamfer curves: from the edge points, geodesics of length d into both faces.
Blend chamfer_blend(const Shell& sh, const Group& g, const NurbsCurve& edge, const NurbsSurface& S1,
                    const NurbsSurface& S2, const std::function<double(double)>& dist, double tol) {
  // interior directions at the midpoint fix the sides (left of each face's loop direction)
  const auto side_of = [&](std::size_t face, const NurbsSurface& S) {
    const ShellEdge& e = sh.edges()[g.edges.front()];
    const Coedge* cu = nullptr;
    for (const auto& u : e.uses)
      if (u.face == face) cu = &u;
    const NurbsCurve& pc = loop_of(sh.faces()[face], cu->loop).curves[cu->index];
    const double tm = 0.5 * (pc.domain_start() + pc.domain_end());
    const auto dd = pc.derivatives(tm, 1);
    const auto D = sh.faces()[face].surface().first_derivatives(dd[0].x, dd[0].y);
    const Vec3 into = cross(cross(D[1], D[2]), dd[1].x * D[1] + dd[1].y * D[2]);
    const CurvePoint cp = closest_point(edge, D[0]);
    const Vec3 T = edge.derivatives(cp.u, 1)[1];
    const SurfacePoint sp = closest_point(S, D[0]);
    const Vec3 dir = cross(S.normal(sp.u, sp.v), T);
    return dot(dir, into) > 0.0 ? 1.0 : -1.0;
  };
  const double sg1 = side_of(g.f1, S1), sg2 = side_of(g.f2, S2);
  // foot points are continued from the nearest sample already computed, so that they stay on
  // one side of the seam of a closed surface (the chain of a closed edge starts on the seam)
  struct Track {
    std::map<double, std::array<double, 2>> feet;
    int steps{8};
  } t1, t2;
  const double s_mid = 0.5 * (edge.domain_start() + edge.domain_end());
  // foot point of the chain point at s on S, marched from the nearest known foot point in small
  // steps (a jump across the seam would send the local search to the wrong side)
  const auto foot = [&](const NurbsSurface& S, Track& tr, double s) {
    if (tr.feet.empty()) {  // seed in the middle of the chain, away from any seam
      const SurfacePoint m0 = closest_point(S, edge.evaluate(s_mid));
      tr.feet[s_mid] = {m0.u, m0.v};
    }
    const double h = (edge.domain_end() - edge.domain_start()) / 128.0;
    SurfacePoint sp{};
    while (true) {
      auto it = tr.feet.lower_bound(s);
      if (it == tr.feet.end() ||
          (it != tr.feet.begin() && s - std::prev(it)->first < it->first - s))
        it = std::prev(it);
      const double gap = s - it->first;
      const double t = std::abs(gap) <= h ? s : it->first + (gap > 0 ? h : -h);
      sp = closest_point_local(S, edge.evaluate(t), it->second[0], it->second[1]);
      tr.feet[t] = {sp.u, sp.v};
      if (t == s) break;
    }
    return sp;
  };
  const auto eval = [&](double s) {
    // the exact edge point (on both surfaces; the chain may be a fitted intersection curve with
    // small kinks at its joints) by alternating projections, and the exact edge direction
    SurfacePoint a = foot(S1, t1, s), b = foot(S2, t2, s);
    for (int it = 0; it < 100 && distance(a.point, b.point) > 1e-14 * (1.0 + norm(a.point)); ++it) {
      a = closest_point_local(S1, b.point, a.u, a.v);
      b = closest_point_local(S2, a.point, b.u, b.v);
    }
    const Vec3 n1 = S1.normal(a.u, a.v), n2 = S2.normal(b.u, b.v);
    Vec3 T = cross(n1, n2);
    if (!(norm(T) > 1e-9)) throw GeometryError("chamfer_edges: the faces are tangent at the edge");
    if (dot(T, edge.derivatives(s, 1)[1]) < 0.0) T = -1.0 * T;
    const double d = dist(s);
    const auto w1 = detail::geodesic_point(S1, a.u, a.v, sg1 * cross(n1, T), d, tol, t1.steps);
    const auto w2 = detail::geodesic_point(S2, b.u, b.v, sg2 * cross(n2, T), d, tol, t2.steps);
    const Vec3 q1{w1[0], w1[1], 0.0}, q2{w2[0], w2[1], 0.0};
    return std::vector<Vec3>{on(S1, q1), on(S2, q2), q1, q2};
  };
  const auto check = [&](const std::vector<NurbsCurve>& c, double s, const std::vector<Vec3>& e) {
    return std::max({distance(c[0].evaluate(s), e[0]), distance(c[1].evaluate(s), e[1]),
                     distance(on(S1, c[2].evaluate(s)), e[0]),
                     distance(on(S2, c[3].evaluate(s)), e[1])});
  };
  // 8 samples per knot span; tiny spans (fitted intersection curves have some) are merged
  std::vector<double> T;
  const auto k = unique_knots(edge.knots()).first;
  const double span_min = 1e-3 * (edge.domain_end() - edge.domain_start());
  for (std::size_t i = 0; i + 1 < k.size(); ++i) {
    if (k[i + 1] - k[i] < span_min) continue;
    for (int j = 0; j < 8; ++j) T.push_back(k[i] + (k[i + 1] - k[i]) * j / 8.0);
  }
  T.push_back(edge.domain_end());
  if (T.front() > edge.domain_start()) T.front() = edge.domain_start();
  auto c = fit_family(T, 4, eval, check, tol, "chamfer_edges");
  NurbsSurface B = ruled(c[0], c[1]).transposed();  // u across (0: face 1), v = edge parameter
  // the curves share the blend's v (ruled() parametrises the rulings' direction on [0, 1])
  const auto dv = B.domain_v();
  for (auto& x : c) x = reparametrize(x, dv[0], dv[1]);
  return {std::move(B), std::move(c[0]), std::move(c[1]), std::move(c[2]), std::move(c[3])};
}

// The 3D chain of the group's edges (one curve, in chain order; closed for a loop).
NurbsCurve edge_chain(const Shell& sh, const Group& g, double tol) {
  std::vector<std::size_t> left = g.edges;
  const auto& E = sh.edges();
  // start at a chain end (or anywhere for a loop)
  std::size_t cur = left.front();
  bool forward = true;
  if (!g.ends.empty()) {
    for (const std::size_t e : left)
      if (E[e].v0 == g.ends.front() || E[e].v1 == g.ends.front()) {
        cur = e;
        forward = E[e].v0 == g.ends.front();
        break;
      }
  }
  // pieces are parametrised proportionally to their length, so that the chain has no jump in
  // its parametric speed at the joints (the blend curves are interpolated in this parameter)
  double len = 0.0;
  const auto piece = [&](const NurbsCurve& c) {
    const double l = arc_length(c);
    NurbsCurve r = reparametrize(c, len, len + l);
    len += l;
    return r;
  };
  NurbsCurve chain = piece(forward ? E[cur].curve : E[cur].curve.reversed());
  std::size_t tip = forward ? E[cur].v1 : E[cur].v0;
  left.erase(std::find(left.begin(), left.end(), cur));
  while (!left.empty()) {
    auto it = std::find_if(left.begin(), left.end(),
                           [&](std::size_t e) { return E[e].v0 == tip || E[e].v1 == tip; });
    if (it == left.end()) throw GeometryError("edge blend: the selected edges do not form a chain");
    const bool fw = E[*it].v0 == tip;
    chain = join(chain, piece(fw ? E[*it].curve : E[*it].curve.reversed()), 10 * tol);
    tip = fw ? E[*it].v1 : E[*it].v0;
    left.erase(it);
  }
  return reparametrize(chain, 0.0, 1.0);
}

// A closed chain rotated to start where it crosses the seam of a closed face surface (the
// boundary loop of such a face is cut there), so that foot points can be followed along it
// without jumping across the seam.
NurbsCurve start_at_seam(const Shell& sh, const Group& g, const NurbsCurve& chain, double tol) {
  for (const std::size_t f : {g.f1, g.f2}) {
    const NurbsSurface& S = sh.faces()[f].surface();
    if (!closed_in(S, true) && !closed_in(S, false)) continue;
    const Coedges c = coedges_in(sh, g, f);
    const TrimLoop& loop = loop_of(sh.faces()[f], c.loop);
    const auto du = S.domain_u(), dv = S.domain_v();
    for (const std::size_t i : c.idx) {
      const Vec3 w = at_start(loop.curves[i]);
      const bool on_seam = std::abs(w.x - du[0]) < 1e-9 || std::abs(w.x - du[1]) < 1e-9 ||
                           std::abs(w.y - dv[0]) < 1e-9 || std::abs(w.y - dv[1]) < 1e-9;
      if (!on_seam) continue;
      const double t = closest_point(chain, on(S, w)).u;
      const double a = chain.domain_start(), b = chain.domain_end();
      if (t - a <= 1e-9 * (b - a) || b - t <= 1e-9 * (b - a)) return chain;  // already there
      auto [p, q] = chain.split(t);
      const NurbsCurve r =
          join(q, reparametrize(p, q.domain_end(), q.domain_end() + (t - a)), 10 * tol);
      return reparametrize(r, 0.0, 1.0);
    }
  }
  return chain;
}

// ------------------------------------------------------------------ curves where a blend ends

// v with Bs(u, v) on E (Newton on the signed distance; the foot point (ue, ve) on E is carried
// along as start value).
double meet_v(const NurbsSurface& Bs, const NurbsSurface& E, double u, double v, double& ue,
              double& ve, const std::string& who) {
  const auto dv = Bs.domain_v();
  for (int it = 0; it < 60; ++it) {
    const auto D = Bs.first_derivatives(u, v);
    const SurfacePoint f = closest_point_local(E, D[0], ue, ve);
    ue = f.u;
    ve = f.v;
    const Vec3 n = E.normal(f.u, f.v);
    const double g = dot(D[0] - f.point, n), dg = dot(D[2], n);
    if (std::abs(g) <= 1e-14 * std::max(1.0, norm(D[0]))) break;
    if (!(std::abs(dg) > 1e-300))
      throw GeometryError(who + ": a blend runs parallel to a surface it ends at");
    v = std::clamp(v - g / dg, dv[0], dv[1]);
  }
  return v;
}

// Exact curve where the blend surface Bs meets the surface E, for u from ua to ub (curve
// parameter w in [0, 1]), starting from v0 at ua: interpolated p-curves on Bs and on E.
struct TrimCurve {
  NurbsCurve on_blend, on_end;
};

// With a hint v1 (the v expected at ub) every point starts from the linear guess between v0 and
// v1 instead of the previous solution: where two branches of the intersection cross at the start
// (two blends touching the same face at one point) this keeps to the branch heading for v1.
TrimCurve trim_curve(const NurbsSurface& Bs, const NurbsSurface& E, double ua, double ub, double v0,
                     double tol, const std::string& who, std::optional<double> v1 = std::nullopt) {
  const SurfacePoint e0 = closest_point(E, Bs.evaluate(ua, v0));
  double v = v0, ue = e0.u, ve = e0.v;
  const auto solve = [&](double w) {
    const double u = ua + (ub - ua) * w;
    if (v1) v = v0 + (*v1 - v0) * w;
    v = meet_v(Bs, E, u, v, ue, ve, who);
    const Vec3 P = Bs.evaluate(u, v);
    const SurfacePoint f = closest_point_local(E, P, ue, ve);
    if (f.distance > 0.1 * tol) {
      throw GeometryError(who + ": a blend does not reach a surface it ends at (u = " + num(u) +
                          ")");
    }
    return std::vector<Vec3>{Vec3{u, v, 0.0}, Vec3{f.u, f.v, 0.0}, P};
  };
  solve(0.0);  // converge the start values at the first end
  const auto eval = [&](double w) { return solve(w); };
  const auto check = [&](const std::vector<NurbsCurve>& c, double w, const std::vector<Vec3>& e) {
    return std::max(distance(on(Bs, c[0].evaluate(w)), e[2]),
                    distance(on(E, c[1].evaluate(w)), e[2]));
  };
  auto c = fit_family(uniform(0.0, 1.0, 8), 2, eval, check, tol, who);
  return {std::move(c[0]), std::move(c[1])};
}

std::vector<NurbsCurve> reversed_chain(const std::vector<NurbsCurve>& pieces) {
  std::vector<NurbsCurve> out;
  for (auto it = pieces.rbegin(); it != pieces.rend(); ++it) out.push_back(it->reversed());
  return out;
}

// ------------------------------------------------------------------ groups and vertices

enum class Kind { Fillet, Chamfer };

// Blend size along a chain (parameter s in [0, 1]): through the given (s, size) points with a
// smooth step between neighbours (quintic, C2: the blend surfaces stay smooth enough for the
// cubic fits), constant outside (periodic for a closed chain).
struct SizeLaw {
  std::vector<std::pair<double, double>> pts;  // sorted by s
  double base{0.0};
  bool closed{false};

  [[nodiscard]] double operator()(double s) const {
    if (pts.empty()) return base;
    if (pts.size() == 1) return pts[0].second;
    const auto step = [](double a, double b, double t) {
      t = std::clamp(t, 0.0, 1.0);
      return a + (b - a) * t * t * t * (10.0 - 15.0 * t + 6.0 * t * t);
    };
    if (closed) {
      s -= std::floor(s);
      const auto& f = pts.front();
      const auto& l = pts.back();
      if (s < f.first || s >= l.first) {  // across the seam
        const double span = f.first + 1.0 - l.first;
        const double t = (s >= l.first ? s - l.first : s + 1.0 - l.first) / span;
        return step(l.second, f.second, t);
      }
    } else {
      if (s <= pts.front().first) return pts.front().second;
      if (s >= pts.back().first) return pts.back().second;
    }
    std::size_t i = 0;
    while (i + 2 < pts.size() && s >= pts[i + 1].first) ++i;
    const double span = pts[i + 1].first - pts[i].first;
    return span <= 0.0 ? pts[i].second
                       : step(pts[i].second, pts[i + 1].second, (s - pts[i].first) / span);
  }
  [[nodiscard]] double max() const {
    double m = pts.empty() ? base : 0.0;
    for (const auto& p : pts) m = std::max(m, p.second);
    return m;
  }
  [[nodiscard]] double min() const {
    double m = pts.empty() ? base : std::numeric_limits<double>::infinity();
    for (const auto& p : pts) m = std::min(m, p.second);
    return m;
  }
};

// A blended group: its blend surface (u = 0 on face f1, u = 1 on face f2, v along the edges),
// and the trimming curves at its end vertices.
struct GroupBlend {
  Group g;
  Blend b;
  std::optional<NurbsCurve> spine;  // fillets: path of the ball centre (parameter v)
  bool convex{true};
  // per end vertex (order of g.ends): connected p-curves on the blend from u = 0 to u = 1
  std::vector<std::vector<NurbsCurve>> trims;
  std::function<double(const Vec3&)> size_at;  // radius / distance at a point of the chain

  [[nodiscard]] bool has(std::size_t f) const { return f == g.f1 || f == g.f2; }
  [[nodiscard]] double u_on(std::size_t f) const { return f == g.f1 ? 0.0 : 1.0; }
  [[nodiscard]] std::size_t other(std::size_t f) const { return f == g.f1 ? g.f2 : g.f1; }
  [[nodiscard]] const NurbsCurve& contact(std::size_t f) const { return f == g.f1 ? b.k1 : b.k2; }
  [[nodiscard]] const NurbsCurve& pcurve(std::size_t f) const { return f == g.f1 ? b.pc1 : b.pc2; }
  // stores the trim at vertex vi given as pieces running from face `from` across the blend
  void set_trim(std::size_t vi, std::vector<NurbsCurve> pieces, std::size_t from) {
    if (u_on(from) != 0.0) pieces = reversed_chain(pieces);
    const auto at = std::find(g.ends.begin(), g.ends.end(), vi);
    trims[static_cast<std::size_t>(at - g.ends.begin())] = std::move(pieces);
  }
};

// What happens to the face loops at a vertex: p-curves inserted between two neighbouring loop
// curves of a face (where a blend ends on it), and points where two contact curves meet.
struct VertexPlan {
  std::map<std::size_t, std::vector<NurbsCurve>> inserts;
  std::map<std::size_t, Vec3> meet;
};

// A new face with a reference direction (sum of the outward normals of its neighbours) that
// fixes its orientation.
struct NewFace {
  TrimmedSurface face;
  Vec3 uv;   // a point inside it (parameter space)
  Vec3 ref;  // its normal there must point this way
};

TrimmedSurface oriented_face(const NewFace& f) {
  const Vec3 n = f.face.surface().normal(f.uv.x, f.uv.y);
  return dot(n, f.ref) < 0.0 ? f.face.flipped() : f.face;
}

const NurbsCurve& edge_pcurve(const Shell& sh, std::size_t e, std::size_t face) {
  for (const auto& u : sh.edges()[e].uses)
    if (u.face == face) return loop_of(sh.faces()[face], u.loop).curves[u.index];
  throw GeometryError("edge blend: internal error (edge not on face)");
}

// Parameter point of vertex V on a face (from the trim curves ending there).
Vec3 uv_at(const Shell& sh, std::size_t face, std::size_t vi) {
  const Vec3 V = sh.vertices()[vi];
  const NurbsSurface& S = sh.faces()[face].surface();
  Vec3 best{};
  double bd = std::numeric_limits<double>::infinity();
  for (std::size_t e = 0; e < sh.edges().size(); ++e) {
    const ShellEdge& E = sh.edges()[e];
    if (E.v0 != vi && E.v1 != vi) continue;
    for (const auto& u : E.uses) {
      if (u.face != face) continue;
      const NurbsCurve& pc = loop_of(sh.faces()[face], u.loop).curves[u.index];
      for (const Vec3& w : {at_start(pc), at_end(pc)})
        if (distance(on(S, w), V) < bd) {
          bd = distance(on(S, w), V);
          best = w;
        }
    }
  }
  if (!(bd < std::numeric_limits<double>::infinity()))
    return Vec3{closest_point(S, V).u, closest_point(S, V).v, 0.0};
  return best;
}

// The crossing of two p-curves (on one face) nearest to `near`.
CurveCurveIntersection crossing(const NurbsCurve& a, const NurbsCurve& b, const Vec3& near,
                                const std::string& who, const std::string& what) {
  // p-curves computed by foot points may end a hair off a seam: retry with a looser tolerance
  auto hits = intersect(a, b, Tolerance{1e-10, tol::kParametric, tol::kAngular});
  if (hits.empty()) hits = intersect(a, b, Tolerance{1e-7, tol::kParametric, tol::kAngular});
  if (hits.empty()) throw GeometryError(who + ": " + what);
  const auto* best = &hits.front();
  for (const auto& h : hits)
    if (distance(h.point, near) < distance(best->point, near)) best = &h;
  return *best;
}

double signed_distance(const NurbsSurface& S, const Vec3& p, double& u, double& v) {
  const SurfacePoint f = closest_point_local(S, p, u, v);
  u = f.u;
  v = f.v;
  return dot(p - f.point, S.normal(f.u, f.v));
}

// Centre of the ball of radius r at signed distance side * r from three surfaces (normals
// outward), by Newton on the three signed distances from the start O.
Vec3 ball_in_corner(const std::array<const NurbsSurface*, 3>& S, double side, double r, Vec3 O,
                    const std::string& who) {
  std::array<SurfacePoint, 3> feet{};
  for (std::size_t q = 0; q < 3; ++q) feet[q] = closest_point(*S[q], O);
  for (int it = 0; it < 50; ++it) {
    std::array<Vec3, 3> n{};
    std::array<double, 3> e{};
    for (std::size_t q = 0; q < 3; ++q) {
      feet[q] = closest_point_local(*S[q], O, feet[q].u, feet[q].v);
      n[q] = S[q]->normal(feet[q].u, feet[q].v);
      e[q] = dot(O - feet[q].point, n[q]) - side * r;
    }
    const double det = dot(n[0], cross(n[1], n[2]));
    if (std::abs(det) <= 1e-12)
      throw GeometryError(who + ": degenerate corner (two faces are tangent)");
    const Vec3 delta =
        (-e[0] * cross(n[1], n[2]) - e[1] * cross(n[2], n[0]) - e[2] * cross(n[0], n[1])) / det;
    O += delta;
    if (norm(delta) <= 1e-15 * std::max(1.0, norm(O))) break;
  }
  return O;
}

class VertexPlanner {
 public:
  VertexPlanner(const Shell& sh, std::vector<GroupBlend>& groups, Kind kind, double size,
                double stol, double btol, std::string who)
      : sh_(sh),
        G_(groups),
        kind_(kind),
        size_(size),
        stol_(stol),
        btol_(btol),
        who_(std::move(who)) {}

  void plan(std::size_t vi, const std::vector<std::size_t>& at) {
    // the faces and edges at the vertex
    std::vector<std::size_t> edges;
    std::set<std::size_t> faces;
    for (std::size_t e = 0; e < sh_.edges().size(); ++e) {
      const ShellEdge& E = sh_.edges()[e];
      if (E.v0 != vi && E.v1 != vi) continue;
      edges.push_back(e);
      for (const auto& u : E.uses) faces.insert(u.face);
    }
    if (faces.size() != 3 || edges.size() != 3)
      throw GeometryError(who_ + ": a blended edge ends at a vertex with " +
                          std::to_string(faces.size()) + " faces and " +
                          std::to_string(edges.size()) +
                          " edges (only vertices of three faces are supported)");
    VertexPlan& P = plans[vi];
    try {
      if (at.size() == 1) {
        single(vi, G_[at[0]], faces, P);
      } else if (at.size() == 2) {
        double_end(vi, G_[at[0]], G_[at[1]], edges, P);
      } else {
        corner(vi, at, faces, P);
      }
    } catch (const GeometryError& e) {
      const Vec3 p = V(vi);
      throw GeometryError(std::string(e.what()) + " [at the vertex (" + num(p.x) + ", " + num(p.y) +
                          ", " + num(p.z) + ")]");
    }
  }

  std::map<std::size_t, VertexPlan> plans;
  std::vector<NewFace> corners;

 private:
  const NurbsSurface& ext(std::size_t f) {
    auto it = ext_.find(f);
    if (it == ext_.end()) it = ext_.emplace(f, extended(sh_.faces()[f].surface(), 3 * size_)).first;
    return it->second;
  }
  Vec3 V(std::size_t vi) const { return sh_.vertices()[vi]; }
  // Faces a and b meet tangentially (same normal) at the vertex.
  bool tangent_at(std::size_t a, std::size_t b, std::size_t vi) const {
    const NurbsSurface& Sa = sh_.faces()[a].surface();
    const NurbsSurface& Sb = sh_.faces()[b].surface();
    const SurfacePoint pa = closest_point(Sa, V(vi)), pb = closest_point(Sb, V(vi));
    const Vec3 na = Sa.normal(pa.u, pa.v), nb = Sb.normal(pb.u, pb.v);
    return norm(cross(na, nb)) <= 1e-6 && dot(na, nb) > 0.0;
  }
  std::size_t edge_between(const std::vector<std::size_t>& edges, std::size_t a, std::size_t b) {
    for (const std::size_t e : edges) {
      const auto& u = sh_.edges()[e].uses;
      if (u.size() == 2 &&
          ((u[0].face == a && u[1].face == b) || (u[0].face == b && u[1].face == a)))
        return e;
    }
    throw GeometryError(who_ + ": no edge between two faces at a blend end");
  }
  // Point where the contact curve of g on face f crosses the edge e (on f).
  Vec3 contact_on_edge(const GroupBlend& g, std::size_t f, std::size_t e, std::size_t vi) {
    const CurveCurveIntersection h = crossing(
        g.pcurve(f), edge_pcurve(sh_, e, f), uv_at(sh_, f, vi), who_,
        "a blend does not reach the neighbouring edge of a face (too large for the face?)");
    return on(sh_.faces()[f].surface(), h.point);
  }

  // One blended edge ends at the vertex: the blend is trimmed by the third face E, whose
  // boundary is rerouted along the end curve.
  void single(std::size_t vi, GroupBlend& g, const std::set<std::size_t>& faces, VertexPlan& P) {
    std::size_t fe = 0;
    for (const std::size_t f : faces)
      if (!g.has(f)) fe = f;
    for (const std::size_t f : {g.g.f1, g.g.f2})
      if (tangent_at(f, fe, vi))
        throw GeometryError(who_ +
                            ": a blended edge continues tangentially into an edge that is not "
                            "blended (select the whole tangent chain)");
    const double v0 = closest_point(g.b.k1, V(vi)).u;
    TrimCurve t = trim_curve(g.b.surface, ext(fe), 0.0, 1.0, v0, btol_, who_);
    P.inserts[fe] = {std::move(t.on_end)};
    g.set_trim(vi, {std::move(t.on_blend)}, g.g.f1);
  }

  // Two blended edges (faces i-j and i-k) meet; the edge j-k is not blended.
  void double_end(std::size_t vi, GroupBlend& A, GroupBlend& B,
                  const std::vector<std::size_t>& edges, VertexPlan& P) {
    const std::size_t i = A.has(B.g.f1) ? B.g.f1 : B.g.f2;
    const std::size_t j = A.other(i), k = B.other(i);
    const std::size_t ejk = edge_between(edges, j, k);
    if (A.convex != B.convex)
      throw GeometryError(who_ +
                          ": a convex and a concave blended edge meet at a vertex "
                          "(not supported)");
    // j and k tangent at the vertex: the blends continue each other (same cross-section)
    double ua = 0.0, va = 0.0, ub = 0.0, vb = 0.0;
    {
      if (tangent_at(j, k, vi)) {
        const Vec3 xa = contact_on_edge(A, j, ejk, vi), xb = contact_on_edge(B, k, ejk, vi);
        va = closest_point(A.contact(j), xa).u;
        vb = closest_point(B.contact(k), xb).u;
        ua = A.u_on(i);
        ub = B.u_on(i);
        const Vec3 pa = A.b.surface.evaluate(ua, va), pb = B.b.surface.evaluate(ub, vb);
        if (distance(pa, pb) > stol_ || distance(xa, xb) > stol_)
          throw GeometryError(who_ + ": the blends of two tangent edges do not meet at a vertex");
        A.set_trim(vi, {make_line({0, va, 0}, {1, va, 0})}, A.g.f1);
        B.set_trim(vi, {make_line({0, vb, 0}, {1, vb, 0})}, B.g.f1);
        P.meet[i] = 0.5 * (pa + pb);
        return;
      }
    }
    // mitre: A and B are trimmed by each other from the point where their contact curves cross
    // on face i; the curve ends on face j or k, and the blend that is left over there is
    // trimmed by that face (it grows into it: an insert in that face's boundary)
    const CurveCurveIntersection h = crossing(A.pcurve(i), B.pcurve(i), uv_at(sh_, i, vi), who_,
                                              "the contact curves of two blends do not cross");
    const Vec3 X = on(sh_.faces()[i].surface(), h.point);
    P.meet[i] = X;
    ua = A.u_on(i);
    const double uj = A.u_on(j);
    va = closest_point(A.contact(i), X).u;
    const NurbsSurface& Sk = ext(k);
    // B extended across its sides as well, so that foot points on it are not clamped to its
    // contact curves (the parameters are kept)
    const NurbsSurface Bx = extended(B.b.surface, size_);
    // march along A ∩ B and watch the side of face k
    double v = va, fu = 0.0, fv = 0.0, ku = 0.0, kv = 0.0;
    {
      const SurfacePoint f0 = closest_point(Bx, X);
      fu = f0.u;
      fv = f0.v;
      const SurfacePoint k0 = closest_point(Sk, X);
      ku = k0.u;
      kv = k0.v;
    }
    const Vec3 qa = contact_on_edge(A, j, ejk, vi);
    const double vq = closest_point(A.contact(j), qa).u;
    const auto at_w = [&](double w, double& vv, double& bu, double& bv, double& su, double& sv) {
      const double u = ua + (uj - ua) * w;
      vv = meet_v(A.b.surface, Bx, u, va + (vq - va) * w, bu, bv, who_);
      return signed_distance(Sk, A.b.surface.evaluate(u, vv), su, sv);
    };
    const double s0 = at_w(0.0, v, fu, fv, ku, kv);
    constexpr int kSteps = 32;
    double w_lo = 0.0, v_lo = v, w_cross = -1.0;
    for (int n = 1; n <= kSteps; ++n) {
      const double w = n / static_cast<double>(kSteps);
      double vv = v_lo, bu = fu, bv = fv, su = ku, sv = kv;
      const double s = at_w(w, vv, bu, bv, su, sv);
      if ((s > 0.0) != (s0 > 0.0)) {
        double a = w_lo, b = w, v_a = v_lo;
        for (int it = 0; it < 60; ++it) {
          const double m = 0.5 * (a + b);
          double vm = v_a, bu2 = bu, bv2 = bv, su2 = su, sv2 = sv;
          const double sm = at_w(m, vm, bu2, bv2, su2, sv2);
          if ((sm > 0.0) == (s0 > 0.0)) {
            a = m;
            v_a = vm;
          } else {
            b = m;
          }
        }
        w_cross = 0.5 * (a + b);
        v_lo = v_a;
        break;
      }
      w_lo = w;
      v_lo = vv;
      fu = bu;
      fv = bv;
      ku = su;
      kv = sv;
    }
    if (w_cross >= 0.0) {
      const double uc = ua + (uj - ua) * w_cross;
      if (distance(A.b.surface.evaluate(uc, v_lo), qa) > stol_) {
        // A runs into face k
        TrimCurve m = trim_curve(A.b.surface, Bx, ua, uc, va, btol_, who_, v_lo);
        const Vec3 e = at_end(m.on_blend);
        TrimCurve p = trim_curve(A.b.surface, Sk, uc, uj, e.y, btol_, who_);
        A.set_trim(vi, {m.on_blend, p.on_blend}, i);
        B.set_trim(vi, {m.on_end}, i);
        P.inserts[k] = {std::move(p.on_end)};
        return;
      }
    }
    TrimCurve m = trim_curve(A.b.surface, Bx, ua, uj, va, btol_, who_, vq);
    const Vec3 M = on(A.b.surface, at_end(m.on_blend));
    const Vec3 qb = contact_on_edge(B, k, ejk, vi);
    if (distance(M, qb) > stol_) {
      // B runs into face j
      const Vec3 e = at_end(m.on_end);
      TrimCurve p = trim_curve(B.b.surface, ext(j), e.x, B.u_on(k), e.y, btol_, who_);
      A.set_trim(vi, {m.on_blend}, i);
      B.set_trim(vi, {m.on_end, p.on_blend}, i);
      P.inserts[j] = {std::move(p.on_end)};
      return;
    }
    A.set_trim(vi, {m.on_blend}, i);
    B.set_trim(vi, {m.on_end}, i);
  }

  // Three blended edges meet: the blends are cut at their cross-sections through the corner,
  // the corner is closed by a spherical triangle (fillets: the ball touching all three faces)
  // or a planar triangle (chamfers).
  void corner(std::size_t vi, const std::vector<std::size_t>& at, const std::set<std::size_t>& fs,
              VertexPlan& P) {
    const std::vector<std::size_t> f(fs.begin(), fs.end());
    for (const std::size_t g : at)
      if (G_[g].convex != G_[at[0]].convex)
        throw GeometryError(who_ +
                            ": convex and concave blended edges meet at a corner "
                            "(not supported)");
    const auto group_of = [&](std::size_t a, std::size_t b) -> GroupBlend& {
      for (const std::size_t g : at)
        if (G_[g].has(a) && G_[g].has(b)) return G_[g];
      throw GeometryError(who_ + ": internal error (corner groups)");
    };
    if (kind_ == Kind::Fillet) {
      fillet_corner(vi, f, group_of, P);
    } else {
      chamfer_corner(vi, f, group_of, P);
    }
  }

  template <class GroupOf>
  void fillet_corner(std::size_t vi, const std::vector<std::size_t>& f, GroupOf& group_of,
                     VertexPlan& P) {
    GroupBlend& A = group_of(f[0], f[1]);
    const double r = A.size_at(V(vi));
    for (std::size_t a = 0; a < 3; ++a)
      for (std::size_t b = a + 1; b < 3; ++b)
        if (std::abs(group_of(f[a], f[b]).size_at(V(vi)) - r) > 1e-9 * std::max(1.0, r))
          throw GeometryError(who_ +
                              ": the fillet radii of three edges meeting at a corner differ");
    const double side = A.convex ? -1.0 : 1.0;  // signed distance of the centre to each face
    // start: the point of A's spine at distance r from the third face, nearest to the vertex
    const NurbsSurface& S2 = sh_.faces()[f[2]].surface();
    const NurbsCurve& sp = *A.spine;
    double su = 0.0, sv = 0.0;
    {
      const SurfacePoint s0 = closest_point(S2, V(vi));
      su = s0.u;
      sv = s0.v;
    }
    const auto h = [&](double t) {
      double a = su, b = sv;
      return signed_distance(S2, sp.evaluate(t), a, b) - side * r;
    };
    const double t0 = sp.domain_start(), t1 = sp.domain_end();
    constexpr int kN = 64;
    double best_t = -1.0, best_d = std::numeric_limits<double>::infinity();
    double prev_t = t0, prev_h = h(t0);
    for (int n = 1; n <= kN; ++n) {
      const double t = t0 + (t1 - t0) * n / kN;
      const double ht = h(t);
      if ((ht > 0.0) != (prev_h > 0.0)) {
        double a = prev_t, b = t, ha = prev_h;
        for (int it = 0; it < 60; ++it) {
          const double m = 0.5 * (a + b), hm = h(m);
          if ((hm > 0.0) == (ha > 0.0)) {
            a = m;
            ha = hm;
          } else {
            b = m;
          }
        }
        const double tc = 0.5 * (a + b);
        const double d = distance(sp.evaluate(tc), V(vi));
        if (d < best_d) {
          best_d = d;
          best_t = tc;
        }
      }
      prev_t = t;
      prev_h = ht;
    }
    if (best_t < 0.0)
      throw GeometryError(who_ + ": no ball of radius " + num(r) + " fits into a corner");
    const Vec3 O =
        ball_in_corner({&ext(f[0]), &ext(f[1]), &ext(f[2])}, side, r, sp.evaluate(best_t), who_);
    std::array<SurfacePoint, 3> feet{};
    for (std::size_t q = 0; q < 3; ++q) feet[q] = closest_point(ext(f[q]), O);
    std::array<Vec3, 3> c{};
    Vec3 ref{};
    for (std::size_t q = 0; q < 3; ++q) {
      feet[q] = closest_point_local(ext(f[q]), O, feet[q].u, feet[q].v);
      c[q] = feet[q].point;
      ref += ext(f[q]).normal(feet[q].u, feet[q].v);
      if (std::abs(feet[q].distance - r) > 0.1 * btol_)
        throw GeometryError(who_ + ": cannot place the corner ball");
      P.meet[f[q]] = c[q];
    }
    // cut the three blends at the cross-sections through the corner ball
    for (std::size_t a = 0; a < 3; ++a)
      for (std::size_t b = a + 1; b < 3; ++b) {
        GroupBlend& g = group_of(f[a], f[b]);
        const double v = closest_point(*g.spine, O).u;
        g.set_trim(vi, {make_line({0, v, 0}, {1, v, 0})}, g.g.f1);
      }
    // spherical triangle c0 c1 c2 on the ball, on a patch of a sphere around its centre
    // direction (poles and seam far away)
    const Vec3 d = normalized((c[0] + c[1] + c[2]) / 3.0 - O);
    Vec3 aux = std::abs(d.x) < 0.9 ? Vec3{1, 0, 0} : Vec3{0, 1, 0};
    const Vec3 ax = normalized(cross(d, aux));
    const Vec3 side_dir = cross(ax, d);
    double lat = 0.0, lon = 0.0;
    for (const Vec3& q : c) {
      const Vec3 w = normalized(q - O);
      lat = std::max(lat, std::abs(std::asin(std::clamp(dot(w, ax), -1.0, 1.0))));
      lon = std::max(lon, std::abs(std::atan2(dot(w, side_dir), dot(w, d))));
    }
    const double margin = 10.0 * std::numbers::pi / 180.0;
    const double beta = lat + margin, alpha = lon + margin;
    if (beta >= 0.5 * std::numbers::pi || alpha >= std::numbers::pi)
      throw GeometryError(who_ + ": corner too wide for a spherical patch");
    const Vec3 e0 = std::cos(alpha) * d - std::sin(alpha) * side_dir;
    const NurbsSurface sphere = revolve(make_arc(O, e0, ax, r, -beta, beta), O, ax, 2.0 * alpha);
    std::vector<NurbsCurve> loop;
    for (std::size_t q = 0; q < 3; ++q) {
      const Vec3 x = normalized(c[q] - O);
      const Vec3 tgt = c[(q + 1) % 3] - O;
      const Vec3 y = normalized(tgt - dot(tgt, x) * x);
      const double th = std::atan2(dot(tgt, y), dot(tgt, x));
      loop.push_back(project_curve(make_arc(O, x, y, r, 0.0, th), sphere, {0, 0, 0}, btol_).pcurve);
    }
    for (std::size_t q = 0; q < 3; ++q)
      loop[q] = snapped(loop[q], at_start(loop[q]), at_start(loop[(q + 1) % 3]));
    const Vec3 w0 = at_start(loop[0]), w1 = at_start(loop[1]), w2 = at_start(loop[2]);
    if (cross(w1 - w0, w2 - w0).z < 0.0) loop = reversed_chain(loop);
    const SurfacePoint mid = closest_point(sphere, O + r * d);
    corners.push_back(
        {TrimmedSurface(sphere, TrimLoop{std::move(loop)}), Vec3{mid.u, mid.v, 0.0}, ref});
  }

  template <class GroupOf>
  void chamfer_corner(std::size_t vi, const std::vector<std::size_t>& f, GroupOf& group_of,
                      VertexPlan& P) {
    // corner points: where the two chamfer contact curves on each face cross
    std::array<Vec3, 3> p{};
    Vec3 ref{};
    for (std::size_t q = 0; q < 3; ++q) {
      const GroupBlend& a = group_of(f[q], f[(q + 1) % 3]);
      const GroupBlend& b = group_of(f[q], f[(q + 2) % 3]);
      const CurveCurveIntersection h =
          crossing(a.pcurve(f[q]), b.pcurve(f[q]), uv_at(sh_, f[q], vi), who_,
                   "the chamfer contact curves do not cross at a corner");
      const NurbsSurface& S = sh_.faces()[f[q]].surface();
      p[q] = on(S, h.point);
      ref += S.normal(h.point.x, h.point.y);
      P.meet[f[q]] = p[q];
    }
    for (std::size_t a = 0; a < 3; ++a)
      for (std::size_t b = a + 1; b < 3; ++b) {
        GroupBlend& g = group_of(f[a], f[b]);
        const double va = closest_point(g.contact(f[a]), p[a]).u;
        const double vb = closest_point(g.contact(f[b]), p[b]).u;
        const double ua = g.u_on(f[a]), ub = g.u_on(f[b]);
        const Vec3 m = g.b.surface.evaluate(0.5 * (ua + ub), 0.5 * (va + vb));
        if (distance(m, 0.5 * (p[a] + p[b])) > btol_)
          throw GeometryError(who_ + ": chamfer corners are only supported between planar faces (" +
                              num(distance(m, 0.5 * (p[a] + p[b]))) + ")");
        g.set_trim(vi, {make_line({ua, va, 0}, {ub, vb, 0})}, f[a]);
      }
    const NurbsSurface tri = make_bilinear(p[0], p[1], p[2], p[1] + p[2] - p[0]);
    TrimLoop loop{{make_line({0, 0, 0}, {1, 0, 0}), make_line({1, 0, 0}, {0, 1, 0}),
                   make_line({0, 1, 0}, {0, 0, 0})}};
    corners.push_back({TrimmedSurface(tri, std::move(loop)), Vec3{0.25, 0.25, 0.0}, ref});
  }

  const Shell& sh_;
  std::vector<GroupBlend>& G_;
  Kind kind_;
  double size_, stol_, btol_;
  std::string who_;
  std::map<std::size_t, NurbsSurface> ext_;
};

// ------------------------------------------------------------------ face loops

// A blended face loop is rebuilt from entries: original trim curves, and runs of a group's
// edges replaced by the group's contact p-curve (oriented along the loop).
struct Entry {
  explicit Entry(NurbsCurve c) : curve(std::move(c)) {}
  NurbsCurve curve;
  long group{-1};
  std::size_t v_start{0}, v_end{0};  // shell vertices at the start / end along the loop
  Vec3 uv_start, uv_end;             // their parameter points
};

struct Junction {
  double ta{0.0}, tb{0.0};      // end parameter on the incoming, start on the outgoing entry
  std::vector<NurbsCurve> mid;  // inserted p-curves (from the incoming to the outgoing entry)
};

class LoopBuilder {
 public:
  LoopBuilder(const Shell& sh, const std::vector<GroupBlend>& G,
              const std::map<std::size_t, long>& owner,
              const std::map<std::size_t, VertexPlan>& plans, std::string who)
      : sh_(sh), G_(G), owner_(owner), plans_(plans), who_(std::move(who)) {}

  // The rebuilt loop l of face f, or nothing if it is untouched.
  std::optional<TrimLoop> rebuild(std::size_t f, std::size_t l) {
    const TrimmedSurface& F = sh_.faces()[f];
    const TrimLoop& loop = loop_of(F, l);
    const std::size_t n = loop.curves.size();
    // edge and direction of every loop curve
    std::vector<std::size_t> edge(n, sh_.edges().size());
    std::vector<bool> rev(n, false);
    for (std::size_t e = 0; e < sh_.edges().size(); ++e)
      for (const auto& u : sh_.edges()[e].uses)
        if (u.face == f && u.loop == l) {
          edge[u.index] = e;
          rev[u.index] = u.reversed;
        }
    std::vector<long> own(n, -1);
    bool touched = false;
    for (std::size_t i = 0; i < n; ++i) {
      if (edge[i] == sh_.edges().size())
        throw GeometryError(who_ + ": internal error (loop curve without edge)");
      const auto it = owner_.find(edge[i]);
      if (it != owner_.end()) {
        own[i] = it->second;
        touched = true;
      }
      const ShellEdge& E = sh_.edges()[edge[i]];
      for (const std::size_t v : {E.v0, E.v1}) {
        const auto p = plans_.find(v);
        if (p != plans_.end() && p->second.inserts.count(f)) touched = true;
      }
    }
    if (!touched) return std::nullopt;
    const auto vstart = [&](std::size_t i) {
      const ShellEdge& E = sh_.edges()[edge[i]];
      return rev[i] ? E.v1 : E.v0;
    };
    const auto vend = [&](std::size_t i) {
      const ShellEdge& E = sh_.edges()[edge[i]];
      return rev[i] ? E.v0 : E.v1;
    };
    // the whole loop is one group's closed run: the closed contact curve becomes the loop
    if (std::all_of(own.begin(), own.end(), [&](long g) { return g >= 0 && g == own[0]; })) {
      const GroupBlend& g = G_[static_cast<std::size_t>(own[0])];
      if (!g.g.ends.empty()) throw GeometryError(who_ + ": internal error (closed run)");
      NurbsCurve R = oriented(g, f, loop.curves[0], g.pcurve(f), F.surface(), true);
      if (distance(at_start(R), at_end(R)) > 1e-7)
        throw GeometryError(who_ + ": the contact curve of a closed edge is not closed on a face");
      return TrimLoop{{snapped(R, at_start(R), at_start(R))}};
    }
    // entries, starting at a boundary between runs
    std::size_t s0 = 0;
    for (std::size_t i = 0; i < n; ++i)
      if (own[i] < 0 || own[i] != own[(i + n - 1) % n]) {
        s0 = i;
        break;
      }
    std::vector<Entry> entries;
    for (std::size_t k = 0; k < n;) {
      const std::size_t i = (s0 + k) % n;
      Entry e(own[i] < 0
                  ? loop.curves[i]
                  : oriented(G_[static_cast<std::size_t>(own[i])], f, loop.curves[i],
                             G_[static_cast<std::size_t>(own[i])].pcurve(f), F.surface(), false));
      e.group = own[i];
      e.v_start = vstart(i);
      e.uv_start = at_start(loop.curves[i]);
      std::size_t last = i;
      ++k;
      if (own[i] >= 0)
        while (k < n && own[(s0 + k) % n] == own[i]) {
          last = (s0 + k) % n;
          ++k;
        }
      e.v_end = vend(last);
      e.uv_end = at_end(loop.curves[last]);
      entries.push_back(std::move(e));
    }
    const std::size_t m = entries.size();
    std::vector<Junction> J(m);  // J[q]: between entries[q - 1] and entries[q]
    for (std::size_t q = 0; q < m; ++q) J[q] = junction(f, entries[(q + m - 1) % m], entries[q]);
    std::vector<NurbsCurve> out;
    for (std::size_t q = 0; q < m; ++q) {
      const Junction& a = J[q];
      const Junction& b = J[(q + 1) % m];
      const NurbsCurve& c = entries[q].curve;
      if (!(b.ta > a.tb))
        throw GeometryError(who_ +
                            ": a blend is too large for a face (its contact curve runs "
                            "past a neighbouring edge)");
      const bool whole = a.tb == c.domain_start() && b.ta == c.domain_end();
      NurbsCurve piece = whole ? c : c.sub_curve(a.tb, b.ta);
      const Vec3 p0 = out.empty() ? c.evaluate(a.tb) : at_end(out.back());
      out.push_back(snapped(piece, p0, c.evaluate(b.ta)));
      for (const NurbsCurve& x : b.mid) out.push_back(snapped(x, at_end(out.back()), at_end(x)));
    }
    // close the loop exactly
    const NurbsCurve& c0 = entries[0].curve;
    out.back() = snapped(out.back(), at_start(out.back()), c0.evaluate(J[0].tb));
    out.front() = snapped(out.front(), at_end(out.back()), at_end(out.front()));
    return TrimLoop{std::move(out)};
  }

 private:
  // The contact p-curve R reoriented to run along the loop (compared with the loop curve c
  // replaced by it, at its middle or start).
  static NurbsCurve oriented(const GroupBlend& g, std::size_t f, const NurbsCurve& c,
                             const NurbsCurve& R, const NurbsSurface& S, bool at_start_point) {
    const double t = at_start_point ? c.domain_start() : 0.5 * (c.domain_start() + c.domain_end());
    const auto d0 = c.derivatives(t, 1);
    const auto D = S.first_derivatives(d0[0].x, d0[0].y);
    const Vec3 T = d0[1].x * D[1] + d0[1].y * D[2];
    const NurbsCurve& k = g.contact(f);
    const CurvePoint cp = closest_point(k, D[0]);
    return dot(k.derivatives(cp.u, 1)[1], T) < 0.0 ? R.reversed() : R;
  }

  Junction junction(std::size_t f, const Entry& a, const Entry& b) {
    Junction J;
    J.ta = a.curve.domain_end();
    J.tb = b.curve.domain_start();
    const VertexPlan* plan = nullptr;
    if (const auto pit = plans_.find(a.v_end); pit != plans_.end()) plan = &pit->second;
    const std::vector<NurbsCurve>* ins = nullptr;
    if (plan != nullptr)
      if (const auto it = plan->inserts.find(f); it != plan->inserts.end()) ins = &it->second;
    if (a.group < 0 && b.group < 0 && ins == nullptr) return J;
    if (ins != nullptr) {
      std::vector<NurbsCurve> mid = *ins;
      const Vec3 s = at_start(mid.front()), e = at_end(mid.back());
      if (closest_point(a.curve, e).distance < closest_point(a.curve, s).distance)
        mid = reversed_chain(mid);
      const CurvePoint pa = closest_point(a.curve, at_start(mid.front()));
      const CurvePoint pb = closest_point(b.curve, at_end(mid.back()));
      if (pa.distance > 1e-6 || pb.distance > 1e-6)
        throw GeometryError(who_ + ": a blend end does not meet the edges of a face");
      J.ta = pa.u;
      J.tb = pb.u;
      J.mid = std::move(mid);
      return J;
    }
    if (a.group >= 0 && b.group >= 0 && plan != nullptr && plan->meet.count(f) != 0) {
      const NurbsSurface& S = sh_.faces()[f].surface();
      const SurfacePoint sp = closest_point_local(S, plan->meet.at(f), a.uv_end.x, a.uv_end.y);
      const Vec3 w{sp.u, sp.v, 0.0};
      J.ta = closest_point(a.curve, w).u;
      J.tb = closest_point(b.curve, w).u;
      return J;
    }
    const CurveCurveIntersection h = crossing(
        a.curve, b.curve, a.uv_end, who_,
        "a blend does not reach the neighbouring edge of a face (too large for the face?)");
    J.ta = h.s;
    J.tb = h.t;
    return J;
  }

  const Shell& sh_;
  const std::vector<GroupBlend>& G_;
  const std::map<std::size_t, long>& owner_;
  const std::map<std::size_t, VertexPlan>& plans_;
  std::string who_;
};

// ------------------------------------------------------------------ driver

Shell blend_edges(const Shell& sh, const std::vector<std::size_t>& edges, double size,
                  const EdgeBlendOptions& opt, Kind kind) {
  const std::string who = kind == Kind::Fillet ? "fillet_edges" : "chamfer_edges";
  if (!(size > 0.0))
    throw GeometryError(who + (kind == Kind::Fillet ? ": radius" : ": distance") +
                        " must be positive");
  if (edges.empty()) return sh;
  const double stol = opt.tolerance > 0.0 ? opt.tolerance : sh.tolerance();
  const double btol = stol / 20.0;
  const std::string label =
      !opt.label.empty() ? opt.label : (kind == Kind::Fillet ? "fillet" : "chamfer");

  // groups by face pair
  std::vector<Group> groups;
  std::map<std::size_t, long> owner;  // edge -> group
  for (const std::size_t e : edges) {
    if (e >= sh.edges().size()) throw GeometryError(who + ": edge index out of range");
    const ShellEdge& E = sh.edges()[e];
    if (E.uses.size() != 2 || E.uses[0].face == E.uses[1].face)
      throw GeometryError(who + ": edge " + std::to_string(e) +
                          " is not shared by exactly two different faces");
    if (owner.count(e) != 0) continue;
    const std::size_t a = E.uses[0].face, b = E.uses[1].face;
    auto it = std::find_if(groups.begin(), groups.end(), [&](const Group& g) {
      return (g.f1 == a && g.f2 == b) || (g.f1 == b && g.f2 == a);
    });
    if (it == groups.end()) {
      groups.push_back(Group{a, b, {}, {}});
      it = std::prev(groups.end());
    }
    it->edges.push_back(e);
    owner[e] = it - groups.begin();
  }
  std::map<std::size_t, std::vector<std::size_t>> at_vertex;  // end vertex -> groups
  for (std::size_t gi = 0; gi < groups.size(); ++gi) {
    Group& g = groups[gi];
    std::map<std::size_t, int> count;
    for (const std::size_t e : g.edges) {
      ++count[sh.edges()[e].v0];
      ++count[sh.edges()[e].v1];
    }
    for (const auto& [v, c] : count)
      if (c == 1) g.ends.push_back(v);
    if (g.ends.size() != 0 && g.ends.size() != 2)
      throw GeometryError(who + ": the edges between two faces must form one chain");
    for (const std::size_t v : g.ends) at_vertex[v].push_back(gi);
  }
  // a vertex where one chain ends and another passes through
  for (const auto& [v, gs] : at_vertex)
    for (std::size_t gi = 0; gi < groups.size(); ++gi) {
      if (std::find(gs.begin(), gs.end(), gi) != gs.end()) continue;
      for (const std::size_t e : groups[gi].edges)
        if (sh.edges()[e].v0 == v || sh.edges()[e].v1 == v)
          throw GeometryError(who +
                              ": a blended edge ends inside another blended edge chain "
                              "(vertex with more than three faces)");
    }

  // blend surfaces (on the face surfaces extended past the edges' ends; a natural extension of
  // a fitted, curved face can bend away quickly, so shorter extensions are tried if it fails)
  std::vector<GroupBlend> G;
  std::vector<bool> used(opt.sizes.size(), false);
  double smax = size;
  for (const auto& p : opt.sizes) {
    if (!(p.size > 0.0))
      throw GeometryError(who + (kind == Kind::Fillet ? ": radii" : ": distances") +
                          " must be positive");
    smax = std::max(smax, p.size);
  }
  // chains and the size along them (in the arc-length fraction: the chain parameter of e.g. a
  // rational circle is not smooth enough in arc length for the cubic fits of the blend)
  struct Prep {
    NurbsCurve chain0;
    bool convex;
    std::shared_ptr<SizeLaw> law;
    std::function<double(double)> frac;
  };
  std::vector<Prep> prep;
  for (const Group& g : groups) {
    NurbsCurve chain0 = edge_chain(sh, g, stol);
    if (g.ends.empty()) chain0 = start_at_seam(sh, g, chain0, stol);
    const auto chain_ptr = std::make_shared<NurbsCurve>(chain0);
    const double total = arc_length(chain0);
    const std::function<double(double)> frac = [chain_ptr, total](double u) {
      const double a = chain_ptr->domain_start(), b = chain_ptr->domain_end();
      return arc_length(*chain_ptr, a, std::clamp(u, a, b)) / total;
    };
    auto law = std::make_shared<SizeLaw>();
    law->base = size;
    law->closed = g.ends.empty();
    for (std::size_t i = 0; i < opt.sizes.size(); ++i) {
      const CurvePoint cp = closest_point(chain0, opt.sizes[i].point);
      if (cp.distance > 10 * stol) continue;
      law->pts.emplace_back(frac(cp.u), opt.sizes[i].size);
      used[i] = true;
    }
    std::sort(law->pts.begin(), law->pts.end());
    prep.push_back({std::move(chain0), is_convex(sh, sh.edges()[g.edges.front()]), law, frac});
  }
  for (std::size_t i = 0; i < used.size(); ++i)
    if (!used[i])
      throw GeometryError(who + ": size point " + std::to_string(i) +
                          " does not lie on a selected edge");
  // corners of three fillets: the ball touching the three faces must fit all three fillets, so
  // a radius given at the corner vertex holds up to the cross-section through that ball
  if (kind == Kind::Fillet && !opt.sizes.empty())
    for (const auto& [v, gs] : at_vertex) {
      if (gs.size() != 3) continue;
      const Vec3 V = sh.vertices()[v];
      const auto at_v = [&](std::size_t gi) {
        const Prep& p = prep[gi];
        return p.frac(closest_point(p.chain0, V).u);
      };
      const double r = (*prep[gs[0]].law)(at_v(gs[0]));
      std::set<std::size_t> fs;
      for (const std::size_t gi : gs) {
        fs.insert(groups[gi].f1);
        fs.insert(groups[gi].f2);
      }
      if (fs.size() != 3) continue;  // reported by the vertex planner
      const std::vector<std::size_t> f(fs.begin(), fs.end());
      Vec3 nsum{};
      std::array<const NurbsSurface*, 3> S{};
      for (std::size_t q = 0; q < 3; ++q) {
        S[q] = &sh.faces()[f[q]].surface();
        const SurfacePoint sp = closest_point(*S[q], V);
        nsum += S[q]->normal(sp.u, sp.v);
      }
      const double side = prep[gs[0]].convex ? -1.0 : 1.0;
      const Vec3 O = ball_in_corner(S, side, r, V + side * r * nsum, who);
      for (const std::size_t gi : gs) {
        SizeLaw& law = *prep[gi].law;
        const double fv = at_v(gi), fo = prep[gi].frac(closest_point(prep[gi].chain0, O).u);
        for (auto& p : law.pts)
          if (std::abs(p.first - fv) <= 1e-9) p.first = fo;
        std::sort(law.pts.begin(), law.pts.end());
        if (std::abs(law(fo) - r) > 1e-9 * std::max(1.0, r))
          throw GeometryError(who +
                              ": the fillet radii of three edges meeting at a corner differ "
                              "(give the same radius at the corner vertex for all three)");
      }
    }

  for (std::size_t gi = 0; gi < groups.size(); ++gi) {
    const Group& g = groups[gi];
    const NurbsCurve& chain0 = prep[gi].chain0;
    const bool convex = prep[gi].convex;
    const auto law = prep[gi].law;
    const auto frac = prep[gi].frac;
    const auto chain_ptr = std::make_shared<NurbsCurve>(chain0);
    const auto size_at = [law, chain_ptr, frac](const Vec3& p) {
      return (*law)(frac(closest_point(*chain_ptr, p).u));
    };
    const double r0 = 0.5 * (law->min() + law->max());  // start values of a variable fillet
    std::optional<GroupBlend> gb;
    for (const double m : {3.0 * smax, 1.5 * smax, 0.75 * smax}) {
      try {
        const NurbsSurface S1 = extended(sh.faces()[g.f1].surface(), m);
        const NurbsSurface S2 = extended(sh.faces()[g.f2].surface(), m);
        std::vector<std::vector<NurbsCurve>> trims(g.ends.size());
        if (kind == Kind::Fillet) {
          SurfaceFilletOptions fo;
          fo.side1 = fo.side2 = convex ? -1 : 1;
          fo.tolerance = btol;
          fo.offset_tolerance = std::max(1e-4 * r0, btol);  // start values only
          fo.near = chain0.evaluate(0.5);
          if (!law->pts.empty()) fo.radius_at = size_at;
          SurfaceFillet f = fillet_surfaces(S1, S2, r0, fo);
          Blend b{std::move(f.surface), std::move(f.contact1), std::move(f.contact2),
                  std::move(f.pcurve1), std::move(f.pcurve2)};
          gb.emplace(
              GroupBlend{g, std::move(b), std::move(f.spine), convex, std::move(trims), size_at});
        } else {
          NurbsCurve chain = chain0;
          if (!g.ends.empty()) {  // run past the end faces
            chain = extended_curve(chain, CurveEnd::Start, m);
            chain = extended_curve(chain, CurveEnd::End, m);
          }
          const auto dist = [law, frac](double s) { return (*law)(frac(s)); };
          gb.emplace(GroupBlend{g, chamfer_blend(sh, g, chain, S1, S2, dist, btol), std::nullopt,
                                convex, std::move(trims), size_at});
        }
        break;
      } catch (const GeometryError&) {
        if (m == 0.75 * smax) throw;
      }
    }
    G.push_back(std::move(*gb));
  }

  // ends
  VertexPlanner planner(sh, G, kind, smax, stol, btol, who);
  for (const auto& [v, gs] : at_vertex) planner.plan(v, gs);

  // faces along the blends and at their ends
  LoopBuilder lb(sh, G, owner, planner.plans, who);
  std::vector<TrimmedSurface> faces;
  for (std::size_t f = 0; f < sh.faces().size(); ++f) {
    const TrimmedSurface& F = sh.faces()[f];
    bool changed = false;
    std::vector<TrimLoop> loops{F.outer()};
    for (const auto& h : F.holes()) loops.push_back(h);
    for (std::size_t l = 0; l < loops.size(); ++l)
      if (auto r = lb.rebuild(f, l)) {
        loops[l] = std::move(*r);
        changed = true;
      }
    if (!changed) {
      faces.push_back(F);
      continue;
    }
    TrimLoop outer = loops.front();
    std::vector<TrimLoop> holes(loops.begin() + 1, loops.end());
    faces.push_back(
        TrimmedSurface(F.surface(), std::move(outer), std::move(holes)).with_label(F.label()));
  }

  // blend faces, trimmed at their ends; oriented like their neighbours
  for (const GroupBlend& gb : G) {
    const NurbsSurface& Bs = gb.b.surface;
    double vref = 0.5 * (Bs.domain_v()[0] + Bs.domain_v()[1]);
    TrimmedSurface face(Bs);
    if (!gb.g.ends.empty()) {
      auto lo = gb.trims[0], hi = gb.trims[1];
      if (at_start(lo.front()).y > at_start(hi.front()).y) std::swap(lo, hi);
      vref = 0.5 * (at_start(lo.front()).y + at_start(hi.front()).y);
      std::vector<NurbsCurve> loop = lo;
      loop.push_back(make_line(at_end(lo.back()), at_end(hi.back())));
      for (const NurbsCurve& c : reversed_chain(hi)) loop.push_back(c);
      loop.push_back(make_line(at_start(hi.front()), at_start(lo.front())));
      for (std::size_t q = 0; q < loop.size(); ++q)
        loop[q] =
            snapped(loop[q], q == 0 ? at_end(loop.back()) : at_end(loop[q - 1]), at_end(loop[q]));
      loop.front() = snapped(loop.front(), at_end(loop.back()), at_end(loop.front()));
      face = TrimmedSurface(Bs, TrimLoop{std::move(loop)});
    }
    const Vec3 q1 = gb.b.pc1.evaluate(vref), q2 = gb.b.pc2.evaluate(vref);
    const Vec3 ref = sh.faces()[gb.g.f1].surface().normal(q1.x, q1.y) +
                     sh.faces()[gb.g.f2].surface().normal(q2.x, q2.y);
    faces.push_back(oriented_face({face, Vec3{0.5, vref, 0.0}, ref}).with_label(label));
  }
  for (const NewFace& c : planner.corners) faces.push_back(oriented_face(c).with_label(label));

  Shell next = Shell::sew(std::move(faces), stol);
  if (!next.is_closed()) {
    const NurbsCurve& c = next.edges()[next.free_edges().front()].curve;
    const Vec3 p = c.evaluate(0.5 * (c.domain_start() + c.domain_end()));
    throw GeometryError(who + ": the result is not closed (" +
                        std::to_string(next.free_edges().size()) + " free edges, e.g. at (" +
                        num(p.x) + ", " + num(p.y) + ", " + num(p.z) + "))");
  }
  if (!next.is_consistently_oriented()) next = next.oriented();
  return next;
}

}  // namespace

std::vector<std::size_t> edges_between(const Shell& shell, const std::string& label_a,
                                       const std::string& label_b) {
  std::vector<std::size_t> out;
  for (std::size_t i = 0; i < shell.edges().size(); ++i) {
    const auto& u = shell.edges()[i].uses;
    if (u.size() != 2) continue;
    const std::string& la = shell.faces()[u[0].face].label();
    const std::string& lb = shell.faces()[u[1].face].label();
    if ((la == label_a && lb == label_b) || (la == label_b && lb == label_a)) out.push_back(i);
  }
  return out;
}

Shell fillet_edges(const Shell& shell, const std::vector<std::size_t>& edges, double radius,
                   const EdgeBlendOptions& options) {
  return blend_edges(shell, edges, radius, options, Kind::Fillet);
}

Shell chamfer_edges(const Shell& shell, const std::vector<std::size_t>& edges, double distance,
                    const EdgeBlendOptions& options) {
  return blend_edges(shell, edges, distance, options, Kind::Chamfer);
}

}  // namespace nf
