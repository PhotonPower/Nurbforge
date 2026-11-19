#include "nurbforge/topology/edge_blend.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <iterator>
#include <limits>
#include <map>
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

// The surface extended by `m` across every side that is not a seam (parameters kept).
NurbsSurface extended(const NurbsSurface& s, double m) {
  NurbsSurface r = s;
  const auto grow = [&](SurfaceSide side) {
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
                    const NurbsSurface& S2, double d, double tol) {
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

// ------------------------------------------------------------------ loop surgery

// Replaces the run of curves `idx` (cyclically contiguous) of the loop by the part of R
// between its crossings with the neighbouring curves; R is reoriented to run along the loop.
// For a run covering the whole loop R itself (closed) becomes the loop.
TrimLoop reroute(const TrimLoop& loop, const std::vector<std::size_t>& idx, NurbsCurve R,
                 const NurbsSurface& S, const NurbsCurve& k3d, const std::string& who) {
  const std::size_t n = loop.curves.size();
  std::set<std::size_t> in(idx.begin(), idx.end());
  // orientation: compare the loop direction with R's parameter via the 3D contact curve
  std::size_t first = n, last = n;
  for (std::size_t i = 0; i < n; ++i)
    if (in.count(i) && !in.count((i + n - 1) % n)) first = i;
  for (std::size_t i = 0; i < n; ++i)
    if (in.count(i) && !in.count((i + 1) % n)) last = i;
  if (in.size() == n) {  // the whole loop
    const Vec3 a = on(S, at_start(loop.curves[0]));
    const auto d0 = loop.curves[0].derivatives(loop.curves[0].domain_start(), 1);
    const auto D = S.first_derivatives(d0[0].x, d0[0].y);
    const Vec3 T = d0[1].x * D[1] + d0[1].y * D[2];
    const CurvePoint cp = closest_point(k3d, a);
    if (dot(k3d.derivatives(cp.u, 1)[1], T) < 0.0) R = R.reversed();
    if (distance(at_start(R), at_end(R)) > 1e-7)
      throw GeometryError(who + ": the contact curve of a closed edge is not closed on a face");
    return TrimLoop{{snapped(R, at_start(R), at_start(R))}};
  }
  if (first == n || last == n)
    throw GeometryError(who + ": the edges are not contiguous in a loop");
  const NurbsCurve& prev = loop.curves[(first + n - 1) % n];
  const NurbsCurve& next = loop.curves[(last + 1) % n];
  if ((first + n - 1) % n == (last + 1) % n)
    throw GeometryError(who + ": a face bounded by the edge and only one other curve");
  // orient R along the loop: compare the loop tangent with the contact curve's tangent in the
  // middle of the run's first curve (its ends may both lie on a seam point)
  {
    const NurbsCurve& c = loop.curves[first];
    const auto d0 = c.derivatives(0.5 * (c.domain_start() + c.domain_end()), 1);
    const auto D = S.first_derivatives(d0[0].x, d0[0].y);
    const Vec3 T = d0[1].x * D[1] + d0[1].y * D[2];
    const CurvePoint cp = closest_point(k3d, D[0]);
    if (dot(k3d.derivatives(cp.u, 1)[1], T) < 0.0) R = R.reversed();
  }
  const Tolerance itol{1e-10, tol::kParametric, tol::kAngular};
  const auto pick = [&](const NurbsCurve& c, const Vec3& near) {
    const auto hits = intersect(R, c, itol);
    if (hits.empty())
      throw GeometryError(who +
                          ": the blend does not reach the neighbouring edge of a face "
                          "(too large for the face?)");
    const auto* best = &hits.front();
    for (const auto& h : hits)
      if (distance(h.point, near) < distance(best->point, near)) best = &h;
    return *best;
  };
  const CurveCurveIntersection hp = pick(prev, at_end(prev));
  const CurveCurveIntersection hn = pick(next, at_start(next));
  if (!(hn.s > hp.s)) throw GeometryError(who + ": blend contact curve runs backwards on a face");
  const Vec3 p0 = prev.evaluate(hp.t), p1 = next.evaluate(hn.t);
  std::vector<NurbsCurve> out;
  for (std::size_t k = 0; k < n; ++k) {
    const std::size_t i = (last + 1 + k) % n;  // start with `next`
    if (in.count(i)) continue;
    if (i == (last + 1) % n) {
      out.push_back(next.sub_curve(hn.t, next.domain_end()));
      out.front() = snapped(out.front(), p1, at_end(next));
    } else if (i == (first + n - 1) % n) {
      NurbsCurve c = prev.sub_curve(prev.domain_start(), hp.t);
      out.push_back(snapped(c, at_start(prev), p0));
    } else {
      out.push_back(loop.curves[i]);
    }
  }
  out.push_back(snapped(R.sub_curve(hp.s, hn.s), p0, p1));
  return TrimLoop{std::move(out)};
}

// ------------------------------------------------------------------ chain ends

// Exact trimming curve of the blend at an end face: for every u the v with B(u, v) on E
// (Newton on the signed distance), interpolated as p-curves on B (u, v(u)) and on E.
struct EndCurve {
  NurbsCurve on_blend, on_end;
};

EndCurve end_curve(const NurbsSurface& B, const NurbsSurface& E, double v0, const Vec3& V,
                   double tol, const std::string& who) {
  const auto dv = B.domain_v();
  const SurfacePoint e0 = closest_point(E, V);
  double v = v0, ue = e0.u, ve = e0.v;
  const auto solve = [&](double u) {
    for (int it = 0; it < 60; ++it) {
      const auto D = B.first_derivatives(u, v);
      const SurfacePoint f = closest_point_local(E, D[0], ue, ve);
      ue = f.u;
      ve = f.v;
      const Vec3 n = E.normal(f.u, f.v);
      const double g = dot(D[0] - f.point, n), dg = dot(D[2], n);
      if (std::abs(g) <= 1e-14 * std::max(1.0, norm(D[0]))) break;
      if (!(std::abs(dg) > 1e-300))
        throw GeometryError(who + ": the blend runs parallel to an end face");
      v = std::clamp(v - g / dg, dv[0], dv[1]);
    }
    const Vec3 P = B.evaluate(u, v);
    const SurfacePoint f = closest_point_local(E, P, ue, ve);
    if (f.distance > 0.1 * tol)
      throw GeometryError(who + ": the blend does not meet the end face (u = " + num(u) + ")");
    return std::vector<Vec3>{Vec3{u, v, 0.0}, Vec3{f.u, f.v, 0.0}, P};
  };
  solve(0.0);  // converge the start values at the first end
  const auto eval = [&](double u) { return solve(u); };
  const auto check = [&](const std::vector<NurbsCurve>& c, double u, const std::vector<Vec3>& e) {
    return std::max(distance(on(B, c[0].evaluate(u)), e[2]),
                    distance(on(E, c[1].evaluate(u)), e[2]));
  };
  auto c = fit_family(uniform(0.0, 1.0, 8), 2, eval, check, tol, who);
  return {std::move(c[0]), std::move(c[1])};
}

// Reroutes face E's boundary around vertex V along the end curve (p-curve on E): the curves
// meeting at V are cut where the end curve starts and ends.
TrimmedSurface reroute_end(const TrimmedSurface& face, const Vec3& V, NurbsCurve I, double tol,
                           const std::string& who) {
  const NurbsSurface& S = face.surface();
  std::vector<TrimLoop> loops{face.outer()};
  for (const auto& h : face.holes()) loops.push_back(h);
  for (auto& loop : loops) {
    const std::size_t n = loop.curves.size();
    for (std::size_t i = 0; i < n; ++i) {
      const NurbsCurve& cin = loop.curves[i];
      const NurbsCurve& cout = loop.curves[(i + 1) % n];
      if (distance(on(S, at_end(cin)), V) > 10 * tol) continue;
      // orient I from c_in to c_out
      if (closest_point(cin, at_end(I)).distance < closest_point(cin, at_start(I)).distance)
        I = I.reversed();
      const CurvePoint a = closest_point(cin, at_start(I)), b = closest_point(cout, at_end(I));
      if (a.distance > 1e-6 || b.distance > 1e-6)
        throw GeometryError(who + ": the blend end does not meet the edges of the end face");
      std::vector<NurbsCurve> out;
      for (std::size_t k = 0; k < n; ++k) {
        const std::size_t j = (i + 2 + k) % n;  // after c_out ... up to c_in
        if (j == i || j == (i + 1) % n) continue;
        out.push_back(loop.curves[j]);
      }
      out.push_back(snapped(cin.sub_curve(cin.domain_start(), a.u), at_start(cin), a.point));
      out.push_back(snapped(I, a.point, b.point));
      out.push_back(snapped(cout.sub_curve(b.u, cout.domain_end()), b.point, at_end(cout)));
      loop = TrimLoop{std::move(out)};
      TrimLoop outer = loops.front();
      std::vector<TrimLoop> holes(loops.begin() + 1, loops.end());
      return TrimmedSurface(S, std::move(outer), std::move(holes)).with_label(face.label());
    }
  }
  throw GeometryError(who + ": the end vertex is not on the end face's boundary");
}

// ------------------------------------------------------------------ driver

struct Selected {
  std::size_t fa, fb;  // face pair (ordered)
  Vec3 mid;            // edge midpoint (identifies the edge after re-sewing)
};

std::vector<Selected> select(const Shell& sh, const std::vector<std::size_t>& edges,
                             const std::string& who) {
  std::vector<Selected> out;
  for (const std::size_t e : edges) {
    if (e >= sh.edges().size()) throw GeometryError(who + ": edge index out of range");
    const ShellEdge& E = sh.edges()[e];
    if (E.uses.size() != 2 || E.uses[0].face == E.uses[1].face)
      throw GeometryError(who + ": edge " + std::to_string(e) +
                          " is not shared by exactly two different faces");
    const auto& c = E.curve;
    out.push_back({std::min(E.uses[0].face, E.uses[1].face),
                   std::max(E.uses[0].face, E.uses[1].face),
                   c.evaluate(0.5 * (c.domain_start() + c.domain_end()))});
  }
  return out;
}

std::size_t find_edge(const Shell& sh, const Vec3& mid, double tol) {
  std::size_t best = sh.edges().size();
  double bd = std::numeric_limits<double>::infinity();
  for (std::size_t i = 0; i < sh.edges().size(); ++i) {
    const double d = closest_point(sh.edges()[i].curve, mid).distance;
    if (d < bd) {
      bd = d;
      best = i;
    }
  }
  if (best == sh.edges().size() || bd > 100 * tol)
    throw GeometryError("edge blend: a selected edge was lost by an earlier blend");
  return best;
}

enum class Kind { Fillet, Chamfer };

Shell blend_edges(const Shell& shell, const std::vector<std::size_t>& edges, double size,
                  const EdgeBlendOptions& opt, Kind kind) {
  const std::string who = kind == Kind::Fillet ? "fillet_edges" : "chamfer_edges";
  if (!(size > 0.0))
    throw GeometryError(who + (kind == Kind::Fillet ? ": radius" : ": distance") +
                        " must be positive");
  if (edges.empty()) return shell;
  const double stol = opt.tolerance > 0.0 ? opt.tolerance : shell.tolerance();
  const double btol = stol / 20.0;
  const std::string label =
      !opt.label.empty() ? opt.label : (kind == Kind::Fillet ? "fillet" : "chamfer");

  // groups by face pair, on the input shell; vertex checks
  const auto sel = select(shell, edges, who);
  std::vector<std::vector<std::size_t>> groups;  // indices into sel
  for (std::size_t i = 0; i < sel.size(); ++i) {
    auto it = std::find_if(groups.begin(), groups.end(), [&](const auto& g) {
      return sel[g.front()].fa == sel[i].fa && sel[g.front()].fb == sel[i].fb;
    });
    if (it == groups.end())
      groups.push_back({i});
    else
      it->push_back(i);
  }
  {
    std::map<std::size_t, std::size_t> owner;  // vertex -> group
    for (std::size_t g = 0; g < groups.size(); ++g)
      for (const std::size_t i : groups[g])
        for (const std::size_t v : {shell.edges()[edges[i]].v0, shell.edges()[edges[i]].v1}) {
          auto [it, fresh] = owner.emplace(v, g);
          if (!fresh && it->second != g)
            throw GeometryError(who +
                                ": selected edges of different face pairs meet at a vertex "
                                "(vertex blends are not supported yet)");
        }
  }

  Shell sh = shell;
  for (const auto& gi : groups) {
    Group g;
    for (const std::size_t i : gi) g.edges.push_back(find_edge(sh, sel[i].mid, stol));
    const ShellEdge& e0 = sh.edges()[g.edges.front()];
    g.f1 = e0.uses[0].face;
    g.f2 = e0.uses[1].face;
    std::map<std::size_t, int> count;
    for (const std::size_t e : g.edges) {
      ++count[sh.edges()[e].v0];
      ++count[sh.edges()[e].v1];
    }
    for (const auto& [v, c] : count)
      if (c == 1) g.ends.push_back(v);
    if (g.ends.size() != 0 && g.ends.size() != 2)
      throw GeometryError(who + ": the edges between two faces must form one chain");
    // closed single edge: v0 == v1 counts twice, no ends

    const TrimmedSurface& F1 = sh.faces()[g.f1];
    const TrimmedSurface& F2 = sh.faces()[g.f2];
    const bool convex = is_convex(sh, e0);
    const double m = 3.0 * size;
    const NurbsSurface S1 = extended(F1.surface(), m), S2 = extended(F2.surface(), m);
    NurbsCurve chain = edge_chain(sh, g, stol);
    if (g.ends.empty()) chain = start_at_seam(sh, g, chain, stol);

    const Blend B = [&]() -> Blend {
      if (kind == Kind::Fillet) {
        SurfaceFilletOptions fo;
        fo.side1 = fo.side2 = convex ? -1 : 1;
        fo.tolerance = btol;
        fo.offset_tolerance = std::max(1e-4 * size, btol);  // start values only
        fo.near = chain.evaluate(0.5);
        SurfaceFillet f = fillet_surfaces(S1, S2, size, fo);
        return {std::move(f.surface), std::move(f.contact1), std::move(f.contact2),
                std::move(f.pcurve1), std::move(f.pcurve2)};
      }
      if (!g.ends.empty()) {  // run past the end faces
        chain = extend(chain, CurveEnd::Start, m);
        chain = extend(chain, CurveEnd::End, m);
      }
      return chamfer_blend(sh, g, chain, S1, S2, size, btol);
    }();

    // faces along the edge
    const Coedges c1 = coedges_in(sh, g, g.f1), c2 = coedges_in(sh, g, g.f2);
    const auto rebuild = [&](const TrimmedSurface& F, const Coedges& c, const NurbsCurve& pc,
                             const NurbsCurve& k) {
      std::vector<TrimLoop> loops{F.outer()};
      for (const auto& h : F.holes()) loops.push_back(h);
      loops[c.loop] = reroute(loops[c.loop], c.idx, pc, F.surface(), k, who);
      TrimLoop outer = loops.front();
      std::vector<TrimLoop> holes(loops.begin() + 1, loops.end());
      return TrimmedSurface(F.surface(), std::move(outer), std::move(holes)).with_label(F.label());
    };
    std::map<std::size_t, TrimmedSurface> replaced;
    replaced.emplace(g.f1, rebuild(F1, c1, B.pc1, B.k1));
    replaced.emplace(g.f2, rebuild(F2, c2, B.pc2, B.k2));

    // ends: trim the blend by the end faces and reroute their boundaries
    TrimmedSurface blend_face(B.surface);
    if (!g.ends.empty()) {
      std::vector<std::pair<double, NurbsCurve>> cuts;  // (v at u = 0, p-curve on B)
      for (const std::size_t vi : g.ends) {
        const Vec3 V = sh.vertices()[vi];
        std::set<std::size_t> faces;
        for (const auto& e : sh.edges())
          if (e.v0 == vi || e.v1 == vi)
            for (const auto& u : e.uses) faces.insert(u.face);
        faces.erase(g.f1);
        faces.erase(g.f2);
        if (faces.size() != 1)
          throw GeometryError(who + ": an end vertex of the edge touches " +
                              std::to_string(faces.size()) +
                              " further faces (exactly one end face is supported)");
        const std::size_t fe = *faces.begin();
        const double v0 = closest_point(B.k1, V).u;
        const EndCurve ec = end_curve(B.surface, sh.faces()[fe].surface(), v0, V, btol, who);
        replaced.insert_or_assign(fe, reroute_end(sh.faces()[fe], V, ec.on_end, stol, who));
        cuts.emplace_back(ec.on_blend.evaluate(0.0).y, ec.on_blend);
      }
      std::sort(cuts.begin(), cuts.end(),
                [](const auto& a, const auto& b) { return a.first < b.first; });
      const NurbsCurve& lo = cuts[0].second;
      const NurbsCurve& hi = cuts[1].second;
      const Vec3 a0 = at_start(lo), a1 = at_end(lo), b0 = at_start(hi), b1 = at_end(hi);
      TrimLoop loop{{lo, make_line(a1, b1), hi.reversed(), make_line(b0, a0)}};
      blend_face = TrimmedSurface(B.surface, std::move(loop));
    }

    std::vector<TrimmedSurface> faces;
    for (std::size_t i = 0; i < sh.faces().size(); ++i) {
      const auto it = replaced.find(i);
      faces.push_back(it != replaced.end() ? it->second : sh.faces()[i]);
    }
    // all other faces keep their (outward) orientation; the blend face is flipped if it does
    // not match its neighbours
    faces.push_back(blend_face.with_label(label));
    Shell next = Shell::sew(faces, stol);
    if (!next.is_closed())
      throw GeometryError(who + ": the result is not closed (" +
                          std::to_string(next.free_edges().size()) + " free edges)");
    if (!next.is_consistently_oriented()) {
      faces.back() = faces.back().flipped();
      next = Shell::sew(std::move(faces), stol);
      if (!next.is_consistently_oriented()) next = next.oriented();
    }
    sh = std::move(next);
  }
  return sh;
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
