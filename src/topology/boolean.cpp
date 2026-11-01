#include "nurbforge/topology/boolean.hpp"

#include <algorithm>
#include <cmath>
#include <optional>

#include "nurbforge/analysis/closest_point.hpp"
#include "nurbforge/analysis/intersect.hpp"
#include "nurbforge/analysis/surface_intersect.hpp"
#include "nurbforge/core/error.hpp"
#include "nurbforge/fitting/fit.hpp"
#include "nurbforge/geometry/primitives.hpp"

namespace nf {

namespace {

struct Box3 {
  Vec3 lo{1e300, 1e300, 1e300}, hi{-1e300, -1e300, -1e300};
  void add(const Vec3& p) {
    lo = {std::min(lo.x, p.x), std::min(lo.y, p.y), std::min(lo.z, p.z)};
    hi = {std::max(hi.x, p.x), std::max(hi.y, p.y), std::max(hi.z, p.z)};
  }
  [[nodiscard]] bool overlaps(const Box3& o, double eps) const {
    return lo.x <= o.hi.x + eps && o.lo.x <= hi.x + eps && lo.y <= o.hi.y + eps && o.lo.y <= hi.y + eps &&
           lo.z <= o.hi.z + eps && o.lo.z <= hi.z + eps;
  }
};

Box3 face_box(const TrimmedSurface& f) {
  Box3 b;
  for (const auto& p : f.surface().control_points()) b.add(p);
  return b;
}

Box3 shell_box(const Shell& s) {
  Box3 b;
  for (const auto& f : s.faces())
    for (const auto& p : f.surface().control_points()) b.add(p);
  return b;
}

// Parameter points inside the trimmed face. From the midpoints of the outer loop's polyline a
// ray is sent inwards (to the left of the counter-clockwise loop); the midpoint between the
// start and the first crossing with any loop lies inside even for very thin pieces (slivers at
// seams or poles), unlike a fixed inward nudge.
// All candidates are returned, widest spot first; the caller checks them with contains() (the
// exact test near curved loops is too expensive to run on every candidate).
std::vector<std::array<double, 2>> interior_points(const TrimmedSurface& f) {
  const auto poly = f.outer().polyline(8);
  std::vector<std::vector<std::array<double, 2>>> loops{poly};
  for (const auto& h : f.holes()) loops.push_back(h.polyline(8));
  std::vector<std::pair<double, std::array<double, 2>>> found;
  for (std::size_t k = 0; k < poly.size(); ++k) {
    const auto& p = poly[k];
    const auto& q = poly[(k + 1) % poly.size()];
    const double ex = q[0] - p[0], ey = q[1] - p[1], len = std::hypot(ex, ey);
    if (len <= 0.0) continue;
    const double mx = 0.5 * (p[0] + q[0]), my = 0.5 * (p[1] + q[1]);
    const double nx = -ey / len, ny = ex / len;  // inward normal
    // first crossing of the ray m + s n (s > 0) with any loop segment
    double s_hit = 1e300;
    for (const auto& lp : loops)
      for (std::size_t i = 0, j = lp.size() - 1; i < lp.size(); j = i++) {
        if (&lp == &loops[0] && (i == (k + 1) % poly.size()) && j == k) continue;  // own segment
        const double ax = lp[j][0], ay = lp[j][1], bx = lp[i][0] - ax, by = lp[i][1] - ay;
        const double den = nx * by - ny * bx;
        if (std::abs(den) < 1e-300) continue;
        const double s = ((ax - mx) * by - (ay - my) * bx) / den;
        const double t = ((ax - mx) * ny - (ay - my) * nx) / den;
        if (s > 1e-12 * len && t >= 0.0 && t <= 1.0) s_hit = std::min(s_hit, s);
      }
    if (s_hit >= 1e300) continue;
    const std::array<double, 2> c{mx + 0.5 * s_hit * nx, my + 0.5 * s_hit * ny};
    found.emplace_back(s_hit, c);
  }
  std::stable_sort(found.begin(), found.end(), [](const auto& x, const auto& y) { return x.first > y.first; });
  std::vector<std::array<double, 2>> out;
  for (const auto& [w, c] : found) out.push_back(c);
  return out;
}

// Distance of a parameter point from the loops of a face (polyline approximation).
double distance_to_loops(const TrimmedSurface& f, double u, double v) {
  double d = 1e300;
  const auto seg = [&](const std::vector<std::array<double, 2>>& poly) {
    for (std::size_t i = 0, j = poly.size() - 1; i < poly.size(); j = i++) {
      const double ax = poly[j][0], ay = poly[j][1], bx = poly[i][0], by = poly[i][1];
      const double ex = bx - ax, ey = by - ay, l2 = ex * ex + ey * ey;
      const double t = l2 > 0.0 ? std::clamp(((u - ax) * ex + (v - ay) * ey) / l2, 0.0, 1.0) : 0.0;
      d = std::min(d, std::hypot(u - ax - t * ex, v - ay - t * ey));
    }
  };
  seg(f.outer().polyline(16));
  for (const auto& h : f.holes()) seg(h.polyline(16));
  return d;
}

// Parameter intervals of an intersection curve on which both p-curves lie inside their faces.
// Where the curve touches a face boundary tangentially (e.g. a circle touching a box edge), it
// runs within the fit accuracy of the boundary over a short stretch and the in/out test flips
// arbitrarily there: short runs that stay within 10 tolerances of either face's boundary are
// taken as touching – inside runs are dropped, gaps between inside runs are closed.
std::vector<std::pair<double, double>> shared_intervals(const IntersectionCurve& c, const TrimmedSurface& fa,
                                                        const TrimmedSurface& fb, double tolerance) {
  const auto inside = [&](double t) {
    const Vec3 w1 = c.pcurve1.evaluate(t), w2 = c.pcurve2.evaluate(t);
    return fa.contains(w1.x, w1.y) && fb.contains(w2.x, w2.y);
  };
  const double t0 = c.curve.domain_start(), t1 = c.curve.domain_end();
  const int n = 256;
  std::vector<std::pair<double, double>> out;
  bool prev = inside(t0);
  double start = t0, tp = t0;
  const auto edge = [&](double lo, double hi, bool lo_in) {
    for (int it = 0; it < 50; ++it) {
      const double m = 0.5 * (lo + hi);
      (inside(m) == lo_in ? lo : hi) = m;
    }
    return 0.5 * (lo + hi);
  };
  for (int i = 1; i <= n; ++i) {
    const double t = t0 + (t1 - t0) * i / n;
    const bool cur = inside(t);
    if (cur != prev) {
      const double e = edge(tp, t, prev);
      if (prev) out.emplace_back(start, e);
      else start = e;
      prev = cur;
    }
    tp = t;
  }
  if (prev) out.emplace_back(start, t1);
  const double tiny = 1e-9 * (t1 - t0);
  out.erase(std::remove_if(out.begin(), out.end(), [&](const auto& iv) { return iv.second - iv.first <= tiny; }),
            out.end());
  // touching runs (short, and within 10 tolerances of a face boundary at their samples)
  std::vector<NurbsCurve> bounds;
  const auto touching = [&](double lo, double hi) {
    if (hi - lo > 0.05 * (t1 - t0)) return false;
    if (bounds.empty())
      for (const TrimmedSurface* f : {&fa, &fb})
        for (auto& loop : f->boundary_curves(0.1 * tolerance))
          for (auto& b : loop) bounds.push_back(std::move(b));
    for (double s : {0.25, 0.5, 0.75}) {
      const Vec3 P = c.curve.evaluate(lo + s * (hi - lo));
      double d = 1e300;
      for (const auto& b : bounds) d = std::min(d, closest_point(b, P).distance);
      if (d > 10.0 * tolerance) return false;
    }
    return true;
  };
  const bool closed = distance(c.curve.evaluate(t0), c.curve.evaluate(t1)) <= tolerance;
  std::vector<std::pair<double, double>> merged;
  for (const auto& iv : out) {
    if (!merged.empty() && touching(merged.back().second, iv.first)) merged.back().second = iv.second;
    else merged.push_back(iv);
  }
  // on a closed curve the gap may run across the curve's start (checked in its two halves)
  if (closed && !merged.empty()) {
    const double gap_lo = merged.back().second, gap_hi = merged.front().first;
    const double len = (t1 - gap_lo) + (gap_hi - t0);
    if (len > tiny && len <= 0.05 * (t1 - t0) && (gap_lo >= t1 - tiny || touching(gap_lo, t1)) &&
        (gap_hi <= t0 + tiny || touching(t0, gap_hi))) {
      merged.back().second = t1;
      merged.front().first = t0;
    }
  }
  std::vector<std::pair<double, double>> kept;
  for (const auto& iv : merged)
    if (!touching(iv.first, iv.second)) kept.push_back(iv);
  return kept;
}

// How two surfaces relate, from samples of each projected onto the other. Samples whose foot
// point lies on the boundary of the other surface's domain say nothing (the point is beyond its
// extent); if every other sample lies on the other surface, the surfaces coincide where they
// overlap (Coincident), if there are no such samples at all they only touch along their domain
// boundaries (Disjoint) – both are degenerate for the surface/surface intersection.
enum class Relation { General, Coincident, Disjoint };

Relation relation(const NurbsSurface& a, const NurbsSurface& b, double tol) {
  int interior = 0;
  for (int dir = 0; dir < 2; ++dir) {
    const NurbsSurface& p = dir == 0 ? b : a;  // sampled
    const NurbsSurface& q = dir == 0 ? a : b;  // projected onto
    const auto pu = p.domain_u(), pv = p.domain_v(), qu = q.domain_u(), qv = q.domain_v();
    const double eu = 1e-7 * (qu[1] - qu[0]), ev = 1e-7 * (qv[1] - qv[0]);
    for (int i = 0; i <= 8; ++i)
      for (int j = 0; j <= 8; ++j) {
        const Vec3 P = p.evaluate(pu[0] + (pu[1] - pu[0]) * i / 8.0, pv[0] + (pv[1] - pv[0]) * j / 8.0);
        const SurfacePoint f = closest_point(q, P);
        if (f.u <= qu[0] + eu || f.u >= qu[1] - eu || f.v <= qv[0] + ev || f.v >= qv[1] - ev) continue;
        if (f.distance > tol) return Relation::General;
        ++interior;
      }
  }
  return interior > 0 ? Relation::Coincident : Relation::Disjoint;
}

// Parameter-space curves on `onto` of those parts of the loops of face `from` that lie on the
// surface `onto` (coincident faces: the other face's boundary becomes a cut). The 3D loop is
// S_from(pc(t)); runs of samples within tolerance of `onto` are located (end points by
// bisection), and each run is interpolated through its projections at the loop parameter,
// doubling the samples until the midpoint error is below 1e-9 of the parameter domain.
std::vector<NurbsCurve> projected_loops(const TrimmedSurface& from, const NurbsSurface& onto, double tol) {
  const auto ou = onto.domain_u(), ov = onto.domain_v();
  const double psize = std::max(ou[1] - ou[0], ov[1] - ov[0]);
  std::vector<NurbsCurve> out;
  const auto handle = [&](const NurbsCurve& pc) {
    const auto X = [&](double t) {
      const Vec3 w = pc.evaluate(t);
      return from.surface().evaluate(w.x, w.y);
    };
    // within tolerance of `onto`; a foot point clamped to the domain boundary only counts if it is
    // (nearly) exact – otherwise runs would reach a tolerance band beyond the surface's edge,
    // where the projection is constant and the parameter curve has a kink
    const double eu = 1e-12 * (ou[1] - ou[0]), ev = 1e-12 * (ov[1] - ov[0]);
    const auto on = [&](double t) {
      const SurfacePoint f = closest_point(onto, X(t));
      if (f.distance > tol) return false;
      const bool clamped = f.u <= ou[0] + eu || f.u >= ou[1] - eu || f.v <= ov[0] + ev || f.v >= ov[1] - ev;
      return !clamped || f.distance <= 1e-3 * tol;
    };
    const double t0 = pc.domain_start(), t1 = pc.domain_end();
    const std::size_t spans = std::max<std::size_t>(1, pc.knots().size() - 2 * (pc.degree() + 1) + 1);
    const int n = static_cast<int>(std::min<std::size_t>(512, 32 * spans));
    std::vector<std::pair<double, double>> runs;
    bool prev = on(t0);
    double start = t0, tp = t0;
    const auto edge = [&](double lo, double hi, bool lo_on) {
      for (int it = 0; it < 50; ++it) {
        const double m = 0.5 * (lo + hi);
        (on(m) == lo_on ? lo : hi) = m;
      }
      return 0.5 * (lo + hi);
    };
    for (int i = 1; i <= n; ++i) {
      const double t = t0 + (t1 - t0) * i / n;
      const bool cur = on(t);
      if (cur != prev) {
        const double e = edge(tp, t, prev);
        if (prev) runs.emplace_back(start, e);
        else start = e;
        prev = cur;
      }
      tp = t;
    }
    if (prev) runs.emplace_back(start, t1);
    for (const auto& [r0, r1] : runs) {
      if (r1 - r0 <= 1e-9 * (t1 - t0)) continue;
      for (std::size_t N = 5;; N = 2 * N - 1) {
        std::vector<double> ts(N);
        std::vector<Vec3> P(N);
        for (std::size_t k = 0; k < N; ++k) {
          ts[k] = k + 1 == N ? r1 : r0 + (r1 - r0) * static_cast<double>(k) / static_cast<double>(N - 1);
          const SurfacePoint f = closest_point(onto, X(ts[k]));
          P[k] = {f.u, f.v, 0.0};
        }
        const NurbsCurve fit = interpolate_curve(P, 3, ts);
        double err = 0.0;
        for (std::size_t k = 0; k + 1 < N; ++k) {
          const double tm = 0.5 * (ts[k] + ts[k + 1]);
          const SurfacePoint f = closest_point(onto, X(tm));
          err = std::max(err, distance(fit.evaluate(tm), Vec3{f.u, f.v, 0.0}));
        }
        if (err <= 1e-8 * psize || N > 1025) {
          out.push_back(fit);
          break;
        }
      }
    }
  };
  for (const auto& c : from.outer().curves) handle(c);
  for (const auto& h : from.holes())
    for (const auto& c : h.curves) handle(c);
  return out;
}

// Drops cuts that run along the face's loops or along an earlier cut (9 samples within 1e-6 of
// the parameter domain): the same boundary often arrives twice – as the projected loop of a
// coincident face and as an intersection curve with an adjacent face – with slightly different
// fits. Nearly coincident curves would make the curve/curve intersection in the arrangement
// subdivide down to its size limit everywhere.
std::vector<NurbsCurve> prune_cuts(const TrimmedSurface& face, std::vector<NurbsCurve> cuts) {
  const auto du = face.surface().domain_u(), dv = face.surface().domain_v();
  const double eps = 1e-6 * std::max(du[1] - du[0], dv[1] - dv[0]);
  std::vector<const NurbsCurve*> ref;
  for (const auto& c : face.outer().curves) ref.push_back(&c);
  for (const auto& h : face.holes())
    for (const auto& c : h.curves) ref.push_back(&c);
  std::vector<NurbsCurve> kept;
  for (auto& c : cuts) {
    bool redundant = true;
    for (int k = 0; k <= 8 && redundant; ++k) {
      const Vec3 p = c.evaluate(c.domain_start() + (c.domain_end() - c.domain_start()) * k / 8.0);
      bool near = false;
      for (const NurbsCurve* r : ref)
        if (closest_point(*r, p).distance <= eps) {
          near = true;
          break;
        }
      for (const auto& q : kept)
        if (!near && closest_point(q, p).distance <= eps) near = true;
      redundant = near;
    }
    if (!redundant) kept.push_back(std::move(c));
  }
  return kept;
}

// Position of a face piece relative to the other shell.
enum class Where { Out, In, OnSame, OnOpposite };

}  // namespace

bool point_in_shell(const Shell& shell, const Vec3& p, double tolerance) {
  const Box3 b = shell_box(shell);
  const double L = 10.0 * (distance(b.lo, b.hi) + distance(p, 0.5 * (b.lo + b.hi)) + 1.0);
  const Vec3 dirs[] = {normalized(Vec3{0.5377, 0.8013, 0.2622}), normalized(Vec3{-0.3112, 0.4419, -0.8414}),
                       normalized(Vec3{0.7071, -0.6018, 0.3712}), normalized(Vec3{-0.8123, -0.2311, 0.5354}),
                       normalized(Vec3{0.1234, -0.9512, -0.2832})};
  int votes_in = 0, votes = 0;
  for (const auto& d : dirs) {
    const NurbsCurve ray = make_line(p, p + L * d);
    int count = 0;
    bool valid = true;
    for (const auto& f : shell.faces()) {
      for (const auto& h : intersect(ray, f.surface(), {tolerance * 1e-3, 1e-10, 1e-9})) {
        if (!f.contains(h.u, h.v)) continue;
        const auto du = f.surface().domain_u(), dv = f.surface().domain_v();
        const double psize = std::max(du[1] - du[0], dv[1] - dv[0]);
        Vec3 n;
        try {
          n = f.surface().normal(h.u, h.v);
        } catch (const GeometryError&) {
          valid = false;
          break;
        }
        if (distance_to_loops(f, h.u, h.v) < 1e-6 * psize || std::abs(dot(n, d)) < 1e-3) {
          valid = false;  // boundary hit or grazing ray: ambiguous
          break;
        }
        ++count;
      }
      if (!valid) break;
    }
    if (!valid) continue;
    ++votes;
    if (count % 2 == 1) ++votes_in;
    if (votes == 3) break;
  }
  if (votes == 0) throw GeometryError("point_in_shell: no unambiguous ray found");
  return 2 * votes_in > votes;
}

Shell boolean(const Shell& a_in, const Shell& b_in, BooleanOp op, double tolerance) {
  if (!a_in.is_closed() || !b_in.is_closed()) throw GeometryError("boolean: both shells must be closed");
  const Shell A = a_in.oriented(), B = b_in.oriented();
  const auto& FA = A.faces();
  const auto& FB = B.faces();
  std::vector<std::vector<NurbsCurve>> cuts_a(FA.size()), cuts_b(FB.size());
  std::vector<std::vector<std::size_t>> partners_a(FA.size()), partners_b(FB.size());  // coincident faces
  std::vector<std::vector<std::size_t>> touching_a(FA.size()), touching_b(FB.size());  // tangential contact
  SurfaceIntersectionOptions opt;
  opt.tol.length = 1e-3 * tolerance;
  for (std::size_t i = 0; i < FA.size(); ++i)
    for (std::size_t j = 0; j < FB.size(); ++j) {
      if (!face_box(FA[i]).overlaps(face_box(FB[j]), tolerance)) continue;
      const Relation rel = relation(FA[i].surface(), FB[j].surface(), tolerance);
      if (rel == Relation::Disjoint) continue;
      if (rel == Relation::Coincident) {
        // each face is cut along the other's boundary; the common part is classified "on"
        for (auto& c : projected_loops(FB[j], FA[i].surface(), tolerance)) cuts_a[i].push_back(std::move(c));
        for (auto& c : projected_loops(FA[i], FB[j].surface(), tolerance)) cuts_b[j].push_back(std::move(c));
        partners_a[i].push_back(j);
        partners_b[j].push_back(i);
        continue;
      }
      const auto res = intersect(FA[i].surface(), FB[j].surface(), opt);
      if (!res.points.empty()) {  // tangential contact (not cut)
        touching_a[i].push_back(j);
        touching_b[j].push_back(i);
      }
      for (const auto& c : res.curves)
        for (const auto& [t0, t1] : shared_intervals(c, FA[i], FB[j], tolerance)) {
          cuts_a[i].push_back(c.pcurve1.sub_curve(t0, t1));
          cuts_b[j].push_back(c.pcurve2.sub_curve(t0, t1));
        }
    }
  // classification of a piece of `face` against the other shell: "on" if its interior point lies
  // on a coincident face of the other shell (same or opposite outward normal), else ray casting.
  // A piece that only touches the other shell (tangential contact, which is not cut) may have
  // its widest interior point on the contact – ray casting is undecided there, so the first
  // candidate clear of the touching faces is used.
  const auto where = [&](const TrimmedSurface& piece, const std::vector<std::size_t>& partners,
                         const std::vector<std::size_t>& touching, const std::vector<TrimmedSurface>& other_faces,
                         const Shell& other) {
    const auto candidates = interior_points(piece);
    const double clear = 100.0 * tolerance;
    Vec3 best;
    double best_gap = -1.0;
    for (std::size_t k = 0; k < candidates.size(); ++k) {
      const auto& q = candidates[k];
      if (!piece.contains(q[0], q[1])) continue;
      const Vec3 P = piece.surface().evaluate(q[0], q[1]);
      for (std::size_t j : partners) {
        const TrimmedSurface& g = other_faces[j];
        const SurfacePoint f = closest_point(g.surface(), P);
        if (f.distance <= 10.0 * tolerance && g.contains(f.u, f.v)) {
          const double c = dot(piece.surface().normal(q[0], q[1]), g.surface().normal(f.u, f.v));
          return c > 0.0 ? Where::OnSame : Where::OnOpposite;
        }
      }
      double gap = 1e300;
      for (std::size_t j : touching) {
        if (gap < clear) break;
        const SurfacePoint f = closest_point(other_faces[j].surface(), P);
        if (other_faces[j].contains(f.u, f.v)) gap = std::min(gap, f.distance);
      }
      if (gap >= clear) return point_in_shell(other, P, tolerance) ? Where::In : Where::Out;
      if (gap > best_gap) {  // touching everywhere: take the point farthest from the other shell
        best_gap = gap;
        best = P;
      }
    }
    if (best_gap < 0.0) throw GeometryError("boolean: cannot find an interior point of a face piece");
    return point_in_shell(other, best, tolerance) ? Where::In : Where::Out;
  };
  std::vector<TrimmedSurface> result;
  for (std::size_t i = 0; i < FA.size(); ++i)
    for (const auto& p : split(FA[i], prune_cuts(FA[i], std::move(cuts_a[i])))) {
      const Where w = where(p, partners_a[i], touching_a[i], FB, B);
      const bool keep = op == BooleanOp::Union          ? (w == Where::Out || w == Where::OnSame)
                        : op == BooleanOp::Intersection ? (w == Where::In || w == Where::OnSame)
                                                        : (w == Where::Out || w == Where::OnOpposite);
      if (keep) result.push_back(p);
    }
  // B pieces on A are never taken: a common face with the same orientation is already taken from
  // A (union, intersection), all other "on" cases drop out
  for (std::size_t j = 0; j < FB.size(); ++j)
    for (const auto& p : split(FB[j], prune_cuts(FB[j], std::move(cuts_b[j])))) {
      const Where w = where(p, partners_b[j], touching_b[j], FA, A);
      if (op == BooleanOp::Union && w == Where::Out) result.push_back(p);
      if (op == BooleanOp::Intersection && w == Where::In) result.push_back(p);
      if (op == BooleanOp::Difference && w == Where::In) result.push_back(p.flipped());
    }
  if (result.empty()) return Shell{};
  // the pieces keep the orientation of their (outward oriented) shells; difference pieces of B
  // are flipped – re-orienting here would turn inner voids outwards
  return Shell::sew(std::move(result), tolerance);
}

}  // namespace nf
