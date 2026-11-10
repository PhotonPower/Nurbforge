#include "nurbforge/analysis/hatch.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <limits>
#include <numbers>
#include <optional>
#include <string>

#include "nurbforge/analysis/closest_point.hpp"
#include "nurbforge/analysis/intersect.hpp"
#include "nurbforge/core/error.hpp"
#include "nurbforge/geometry/primitives.hpp"
#include "nurbforge/modeling/offset.hpp"
#include "nurbforge/nurbs/knot_vector.hpp"

namespace nf {

namespace {

// ------------------------------------------------------------------ helpers

double cross2(const Vec3& a, const Vec3& b) { return a.x * b.y - a.y * b.x; }
double dot2(const Vec3& a, const Vec3& b) { return a.x * b.x + a.y * b.y; }

std::string fmt(double v) {
  char buf[32];
  std::snprintf(buf, sizeof buf, "%.6g", v);
  return buf;
}

/// Signed area of a loop in the xy plane (Green, Gauss-Legendre per span).
double loop_area(const std::vector<NurbsCurve>& curves) {
  static const double gx[8] = {-0.9602898564975363, -0.7966664774136267, -0.5255324099163290,
                               -0.1834346424956498, 0.1834346424956498,  0.5255324099163290,
                               0.7966664774136267,  0.9602898564975363};
  static const double gw[8] = {0.1012285362903763, 0.2223810344533745, 0.3137066458778873,
                               0.3626837833783620, 0.3626837833783620, 0.3137066458778873,
                               0.2223810344533745, 0.1012285362903763};
  double a = 0.0;
  for (const auto& c : curves) {
    const auto k = unique_knots(c.knots()).first;
    for (std::size_t i = 0; i + 1 < k.size(); ++i) {
      if (k[i + 1] <= k[i] || k[i + 1] <= c.domain_start() || k[i] >= c.domain_end()) continue;
      for (int q = 0; q < 4; ++q) {
        const double a0 = k[i] + (k[i + 1] - k[i]) * q / 4.0,
                     a1 = k[i] + (k[i + 1] - k[i]) * (q + 1) / 4.0;
        for (int g = 0; g < 8; ++g) {
          const double t = 0.5 * (a0 + a1) + 0.5 * (a1 - a0) * gx[g];
          const auto d = c.derivatives(t, 1);
          a += 0.25 * (a1 - a0) * gw[g] * (d[0].x * d[1].y - d[0].y * d[1].x);
        }
      }
    }
  }
  return a;
}

// ------------------------------------------------------------------ line/loop crossings

/// Bezier piece of a loop curve: heights h_i = n . P_i of the control points and the weights.
struct Bez {
  std::size_t curve{0};
  double u0{0.0}, u1{0.0};  ///< parameter interval on the loop curve
  std::vector<double> h, w;
  double hmin{0.0}, hmax{0.0};
};

/// Run of constant state along a loop: -1 below the line (n . p < c - eps), +1 above
/// (n . p > c + eps), 0 within the tolerance band or undecided at the subdivision limit.
struct Run {
  int state{0};
  std::size_t curve{0};
  double u0{0.0}, u1{0.0};
};

void push_run(std::vector<Run>& runs, int state, std::size_t curve, double u0, double u1) {
  if (!runs.empty() && runs.back().state == state && runs.back().curve == curve &&
      runs.back().u1 >= u0) {
    runs.back().u1 = u1;
    return;
  }
  runs.push_back({state, curve, u0, u1});
}

/// de Casteljau split at s of homogeneous values (N_i = w_i a_i, W_i = w_i).
void split_bez(const std::vector<double>& N, const std::vector<double>& W, double s,
               std::vector<double>& NL, std::vector<double>& WL, std::vector<double>& NR,
               std::vector<double>& WR) {
  const std::size_t n = N.size();
  std::vector<double> a = N, b = W;
  NL.resize(n);
  WL.resize(n);
  NR.resize(n);
  WR.resize(n);
  NL[0] = a[0];
  WL[0] = b[0];
  NR[n - 1] = a[n - 1];
  WR[n - 1] = b[n - 1];
  for (std::size_t r = 1; r < n; ++r) {
    for (std::size_t i = 0; i + r < n; ++i) {
      a[i] = (1 - s) * a[i] + s * a[i + 1];
      b[i] = (1 - s) * b[i] + s * b[i + 1];
    }
    NL[r] = a[0];
    WL[r] = b[0];
    NR[n - 1 - r] = a[n - 1 - r];
    WR[n - 1 - r] = b[n - 1 - r];
  }
}

/// Numerator sum N_i B_i(s) (its sign is the sign of the height difference).
double eval_num(const std::vector<double>& N, double s) {
  std::vector<double> a = N;
  for (std::size_t r = 1; r < a.size(); ++r)
    for (std::size_t i = 0; i + r < a.size(); ++i) a[i] = (1 - s) * a[i] + s * a[i + 1];
  return a[0];
}

/// Simple root of the numerator in (0, 1) with N(0) and N(1) of opposite sign (Illinois).
double root_num(const std::vector<double>& N) {
  double a = 0.0, b = 1.0, fa = N.front(), fb = N.back();
  int side = 0;
  for (int it = 0; it < 200 && b - a > 1e-16; ++it) {
    const double x = (a * fb - b * fa) / (fb - fa);
    const double fx = eval_num(N, x);
    if (fx == 0.0) return x;
    if ((fx < 0) == (fa < 0)) {
      a = x;
      fa = fx;
      if (side == -1) fb *= 0.5;
      side = -1;
    } else {
      b = x;
      fb = fx;
      if (side == 1) fa *= 0.5;
      side = 1;
    }
  }
  return 0.5 * (a + b);
}

/// Classifies the Bezier piece with height differences a_i = h_i - c on the local interval
/// [s0, s1] of the piece into runs (recursive subdivision; heights of a rational piece are convex
/// combinations of the control heights).
void classify(const std::vector<double>& N, const std::vector<double>& W, double s0, double s1,
              double eps, int depth, const Bez& bz, std::vector<Run>& runs) {
  const auto u_of = [&](double s) { return bz.u0 + s * (bz.u1 - bz.u0); };
  double amin = std::numeric_limits<double>::max(), amax = -amin;
  int changes = 0;
  int prev = 0;
  for (std::size_t i = 0; i < N.size(); ++i) {
    const double a = N[i] / W[i];
    amin = std::min(amin, a);
    amax = std::max(amax, a);
    const int sg = N[i] > 0 ? 1 : (N[i] < 0 ? -1 : 0);
    if (sg != 0) {
      if (prev != 0 && sg != prev) ++changes;
      prev = sg;
    }
  }
  if (amin > eps) return push_run(runs, 1, bz.curve, u_of(s0), u_of(s1));
  if (amax < -eps) return push_run(runs, -1, bz.curve, u_of(s0), u_of(s1));
  if (amin >= -eps && amax <= eps) return push_run(runs, 0, bz.curve, u_of(s0), u_of(s1));
  const double a0 = N.front() / W.front(), a1 = N.back() / W.back();
  if (changes == 1 && std::abs(a0) > eps && std::abs(a1) > eps) {
    // exactly one (simple) root, strict states on both ends: a transversal crossing
    const double r = s0 + (s1 - s0) * root_num(N);
    push_run(runs, a0 < 0 ? -1 : 1, bz.curve, u_of(s0), u_of(r));
    push_run(runs, a1 < 0 ? -1 : 1, bz.curve, u_of(r), u_of(s1));
    return;
  }
  if (depth >= 64 || s1 - s0 < 1e-15) return push_run(runs, 0, bz.curve, u_of(s0), u_of(s1));
  std::vector<double> NL, WL, NR, WR;
  split_bez(N, W, 0.5, NL, WL, NR, WR);
  const double sm = 0.5 * (s0 + s1);
  classify(NL, WL, s0, sm, eps, depth + 1, bz, runs);
  classify(NR, WR, sm, s1, eps, depth + 1, bz, runs);
}

struct Crossing {
  double x{0.0};  ///< position along the direction
  int wind{0};    ///< +1 entering the region along +d, -1 leaving it
  std::size_t loop{0};
};

/// Crossings of the line n . p = c with one loop (see the header for the contact rules).
void loop_crossings(const std::vector<Bez>& bez, const std::vector<NurbsCurve>& curves,
                    std::size_t loop, double c, double eps, const Vec3& d, std::vector<Run>& runs,
                    std::vector<Crossing>& out) {
  runs.clear();
  std::vector<double> N, W;
  for (const auto& bz : bez) {
    if (bz.hmin - c > eps) {
      push_run(runs, 1, bz.curve, bz.u0, bz.u1);
      continue;
    }
    if (bz.hmax - c < -eps) {
      push_run(runs, -1, bz.curve, bz.u0, bz.u1);
      continue;
    }
    N.resize(bz.h.size());
    W = bz.w;
    for (std::size_t i = 0; i < N.size(); ++i) N[i] = bz.w[i] * (bz.h[i] - c);
    classify(N, W, 0.0, 1.0, eps, 0, bz, runs);
  }
  std::size_t first = runs.size();
  for (std::size_t i = 0; i < runs.size(); ++i)
    if (runs[i].state != 0) {
      first = i;
      break;
    }
  if (first == runs.size()) return;  // the whole loop lies in the band
  std::size_t last = first;
  for (std::size_t k = 1; k <= runs.size(); ++k) {
    const std::size_t j = (first + k) % runs.size();
    if (runs[j].state == 0) continue;
    if (runs[j].state != runs[last].state) {
      // crossing between runs[last] and runs[j]: placed at the end of the band adjacent to the
      // side above the line (exactly at the root for a transversal crossing)
      const bool up = runs[last].state < 0;
      const Run& r = up ? runs[j] : runs[last];
      const Vec3 p = curves[r.curve].evaluate(up ? r.u0 : r.u1);
      out.push_back({dot2(d, p), up ? -1 : 1, loop});
    }
    last = j;
  }
}

// ------------------------------------------------------------------ contours

/// Bezier pieces of all loops with heights along n.
std::vector<std::vector<Bez>> make_bez(const std::vector<std::vector<NurbsCurve>>& loops,
                                       const Vec3& n) {
  std::vector<std::vector<Bez>> bez(loops.size());
  for (std::size_t li = 0; li < loops.size(); ++li)
    for (std::size_t ci = 0; ci < loops[li].size(); ++ci)
      for (const auto& b : loops[li][ci].decompose_bezier()) {
        Bez z;
        z.curve = ci;
        z.u0 = b.domain_start();
        z.u1 = b.domain_end();
        if (!(z.u1 > z.u0)) continue;
        z.w = b.weights();
        z.h.reserve(b.num_control_points());
        for (const auto& p : b.control_points()) z.h.push_back(dot2(n, p));
        z.hmin = *std::min_element(z.h.begin(), z.h.end());
        z.hmax = *std::max_element(z.h.begin(), z.h.end());
        bez[li].push_back(std::move(z));
      }
  return bez;
}

/// Splits a curve at interior knots of full multiplicity (possible tangent discontinuities).
std::vector<NurbsCurve> split_at_c0_knots(const NurbsCurve& c) {
  const auto [k, m] = unique_knots(c.knots());
  std::vector<double> cuts{c.domain_start()};
  for (std::size_t i = 0; i < k.size(); ++i)
    if (k[i] > c.domain_start() && k[i] < c.domain_end() &&
        static_cast<std::size_t>(m[i]) >= c.degree())
      cuts.push_back(k[i]);
  cuts.push_back(c.domain_end());
  if (cuts.size() == 2) return {c};
  std::vector<NurbsCurve> out;
  for (std::size_t i = 0; i + 1 < cuts.size(); ++i)
    out.push_back(c.sub_curve(cuts[i], cuts[i + 1]));
  return out;
}

Vec3 unit_xy(const Vec3& v) {
  const double l = std::hypot(v.x, v.y);
  return l > 0 ? Vec3{v.x / l, v.y / l, 0.0} : Vec3{0, 0, 0};
}

Vec3 start_tangent(const NurbsCurve& c) { return unit_xy(c.derivatives(c.domain_start(), 1)[1]); }
Vec3 end_tangent(const NurbsCurve& c) { return unit_xy(c.derivatives(c.domain_end(), 1)[1]); }

/// Signed curvature in the xy plane (> 0 turning left).
double signed_curvature(const NurbsCurve& c, double u) {
  const auto d = c.derivatives(u, 2);
  const double s = std::hypot(d[1].x, d[1].y);
  return s > 0 ? cross2(d[1], d[2]) / (s * s * s) : 0.0;
}

/// Parameter intervals of c on which d * signed curvature < limit (the offset by d to the left
/// is regular there); the rest would give cusps and lies in the swallowtail of the offset.
std::vector<std::pair<double, double>> regular_intervals(const NurbsCurve& c, double d,
                                                         double limit, int samples) {
  const auto g = [&](double u) { return limit - d * signed_curvature(c, u); };
  const auto k = unique_knots(c.knots()).first;
  std::vector<std::pair<double, double>> out;
  bool open = false;
  double start = 0.0;
  for (std::size_t i = 0; i + 1 < k.size(); ++i) {
    if (k[i + 1] <= k[i] || k[i + 1] <= c.domain_start() || k[i] >= c.domain_end()) continue;
    double ua = k[i], ga = g(ua);
    if (i == 0 || !open) {
      if (ga > 0 && !open) {
        open = true;
        start = ua;
      }
    }
    for (int s = 1; s <= samples; ++s) {
      const double ub = k[i] + (k[i + 1] - k[i]) * s / samples, gb = g(ub);
      if ((ga > 0) != (gb > 0)) {
        double a = ua, b = ub;  // bisection of the sign change
        for (int it = 0; it < 60; ++it) {
          const double m = 0.5 * (a + b);
          if ((g(m) > 0) == (ga > 0))
            a = m;
          else
            b = m;
        }
        const double r = 0.5 * (a + b);
        if (ga > 0) {
          if (r > start) out.emplace_back(start, r);
          open = false;
        } else {
          open = true;
          start = r;
        }
      }
      ua = ub;
      ga = gb;
    }
  }
  if (open && c.domain_end() > start) out.emplace_back(start, c.domain_end());
  return out;
}

/// Splits a curve into pieces whose tangent turns by at most pi/2 (such pieces are simple, so
/// all self-intersections of the offset are found as intersections between pieces).
void split_by_turning(const NurbsCurve& c, std::vector<NurbsCurve>& out) {
  const auto k = unique_knots(c.knots()).first;
  std::vector<double> cuts{c.domain_start()};
  double turned = 0.0;
  Vec3 prev = start_tangent(c);
  for (std::size_t i = 0; i + 1 < k.size(); ++i) {
    if (k[i + 1] <= k[i] || k[i + 1] <= c.domain_start() || k[i] >= c.domain_end()) continue;
    for (int s = 1; s <= 8; ++s) {
      const double u = k[i] + (k[i + 1] - k[i]) * s / 8;
      const Vec3 t = unit_xy(c.derivatives(u, 1)[1]);
      turned += std::abs(std::atan2(cross2(prev, t), dot2(prev, t)));
      prev = t;
      if (turned > 0.5 * std::numbers::pi && u < c.domain_end()) {
        cuts.push_back(u);
        turned = 0.0;
      }
    }
  }
  cuts.push_back(c.domain_end());
  for (std::size_t i = 0; i + 1 < cuts.size(); ++i)
    if (cuts[i + 1] > cuts[i])
      out.push_back(cuts.size() == 2 ? c : c.sub_curve(cuts[i], cuts[i + 1]));
}

struct Raw {
  NurbsCurve curve;
  std::size_t loop{0};
};

/// Raw inward offset of every loop by d: regular parts of the piece offsets plus arcs about the
/// corner at reflex corners; at convex corners and cusps the raw offset overlaps itself.
void raw_offsets(const std::vector<std::vector<NurbsCurve>>& loops, double d, double fit_tol,
                 std::vector<Raw>& raw) {
  for (std::size_t li = 0; li < loops.size(); ++li) {
    std::vector<NurbsCurve> pieces;
    for (const auto& c : loops[li])
      for (auto& p : split_at_c0_knots(c)) pieces.push_back(std::move(p));
    const std::size_t n = pieces.size();
    for (std::size_t i = 0; i < n; ++i) {
      const auto& c = pieces[i];
      for (const int samples : {16, 128}) {
        try {
          std::vector<NurbsCurve> offs;
          for (const auto& [a, b] : regular_intervals(c, d, 1.0 - 1e-6, samples)) {
            const auto part = a > c.domain_start() || b < c.domain_end() ? c.sub_curve(a, b) : c;
            offs.push_back(offset_curve(part, d, {0, 0, 1}, fit_tol));
          }
          for (auto& o : offs) {
            std::vector<NurbsCurve> parts;
            split_by_turning(o, parts);
            for (auto& q : parts) raw.push_back({std::move(q), li});
          }
          break;
        } catch (const GeometryError&) {
          if (samples == 128) throw;
        }
      }
      // corner to the next piece
      const std::size_t j = (i + 1) % n;
      const Vec3 ta = end_tangent(pieces[i]), tb = start_tangent(pieces[j]);
      const double cr = cross2(ta, tb), dt = dot2(ta, tb);
      if (cr < -1e-9 || (std::abs(cr) <= 1e-9 && dt < 0)) {
        // reflex corner (right turn): arc about the corner from left(ta) clockwise to left(tb)
        const Vec3 v = pieces[i].evaluate(pieces[i].domain_end());
        const double phi = std::atan2(std::abs(cr), dt);
        std::vector<NurbsCurve> parts;
        split_by_turning(make_arc(v, {-ta.y, ta.x, 0}, {ta.x, ta.y, 0}, d, 0.0, phi), parts);
        for (auto& q : parts) raw.push_back({std::move(q), li});
      }
    }
  }
}

struct Fragment {
  NurbsCurve curve;
  std::size_t loop{0};
  Vec3 a, b, ta, tb;
  double length{0.0};
  bool used{false};
};

double chord_length(const NurbsCurve& c) {
  double l = 0.0;
  Vec3 p = c.evaluate(c.domain_start());
  for (int i = 1; i <= 16; ++i) {
    const Vec3 q = c.evaluate(c.domain_start() + (c.domain_end() - c.domain_start()) * i / 16);
    l += distance(p, q);
    p = q;
  }
  return l;
}

/// Winding number of the loops about p (crossings of the line y = p.y left of p).
int winding(const std::vector<std::vector<Bez>>& bez_y,
            const std::vector<std::vector<NurbsCurve>>& loops, const Vec3& p, double eps) {
  std::vector<Run> runs;
  std::vector<Crossing> xs;
  for (std::size_t li = 0; li < loops.size(); ++li)
    loop_crossings(bez_y[li], loops[li], li, p.y, eps, {1, 0, 0}, runs, xs);
  int w = 0;
  for (const auto& x : xs)
    if (x.x < p.x) w += x.wind;
  return w;
}

/// Inward offset of the region bounded by the loops at distance d: the raw offsets are split at
/// all their intersections, the fragments farther than d from the boundary (and inside the
/// region) are kept and chained into closed contours (at junctions the leftmost continuation).
std::vector<HatchContour> contours_at(const std::vector<std::vector<NurbsCurve>>& loops,
                                      const std::vector<std::vector<Bez>>& bez_y, double d,
                                      double fit_tol, double eps, double join_gap) {
  std::vector<Raw> raw;
  try {
    raw_offsets(loops, d, fit_tol, raw);
  } catch (const GeometryError& e) {
    throw GeometryError("hatch: contour at distance " + fmt(d) + ": " + e.what());
  }
  const double ctol = std::max(100.0 * fit_tol, 10.0 * join_gap);
  const std::size_t n = raw.size();
  std::vector<std::array<double, 4>> box(n);
  std::vector<Vec3> rs(n), re(n);
  for (std::size_t i = 0; i < n; ++i) {
    box[i] = {1e300, 1e300, -1e300, -1e300};
    for (const auto& q : raw[i].curve.control_points())
      box[i] = {std::min(box[i][0], q.x), std::min(box[i][1], q.y), std::max(box[i][2], q.x),
                std::max(box[i][3], q.y)};
    rs[i] = raw[i].curve.evaluate(raw[i].curve.domain_start());
    re[i] = raw[i].curve.evaluate(raw[i].curve.domain_end());
  }
  // split points at all intersections
  std::vector<std::vector<double>> cuts(n);
  for (std::size_t i = 0; i < n; ++i)
    for (std::size_t j = i + 1; j < n; ++j) {
      if (box[i][0] > box[j][2] + ctol || box[j][0] > box[i][2] + ctol ||
          box[i][1] > box[j][3] + ctol || box[j][1] > box[i][3] + ctol)
        continue;
      for (const auto& x : intersect(raw[i].curve, raw[j].curve, Tolerance{fit_tol, 1e-12, 1e-9})) {
        // the end of one piece touching the end of the other: a connection, not a crossing
        const auto joint = [&](const Vec3& p, const Vec3& q) {
          return distance(x.point, p) <= ctol && distance(x.point, q) <= ctol &&
                 distance(p, q) <= ctol;
        };
        if (joint(re[i], rs[j]) || joint(rs[i], re[j]) || joint(re[i], re[j]) ||
            joint(rs[i], rs[j]))
          continue;
        cuts[i].push_back(x.s);  // fragments shorter than ctol next to an end are dropped
        cuts[j].push_back(x.t);
      }
    }
  // fragments; keep those at distance >= d from the boundary and inside the region
  std::vector<std::array<double, 4>> lbox;
  std::vector<const NurbsCurve*> lcurve;
  for (const auto& l : loops)
    for (const auto& c : l) {
      std::array<double, 4> b{1e300, 1e300, -1e300, -1e300};
      for (const auto& q : c.control_points())
        b = {std::min(b[0], q.x), std::min(b[1], q.y), std::max(b[2], q.x), std::max(b[3], q.y)};
      lbox.push_back(b);
      lcurve.push_back(&c);
    }
  const double dmin = d - 2.0 * fit_tol - join_gap;
  std::vector<Fragment> frags;
  for (std::size_t i = 0; i < n; ++i) {
    const auto& c = raw[i].curve;
    auto& cu = cuts[i];
    cu.push_back(c.domain_start());
    cu.push_back(c.domain_end());
    std::sort(cu.begin(), cu.end());
    const double ptol = 1e-12 * (c.domain_end() - c.domain_start());
    for (std::size_t k = 0; k + 1 < cu.size(); ++k) {
      if (cu[k + 1] - cu[k] <= ptol) continue;
      const double um = 0.5 * (cu[k] + cu[k + 1]);
      const Vec3 m = c.evaluate(um);
      bool keep = true;
      for (std::size_t q = 0; q < lcurve.size() && keep; ++q) {
        const auto& b = lbox[q];
        const double bx = std::max({b[0] - m.x, 0.0, m.x - b[2]}),
                     by = std::max({b[1] - m.y, 0.0, m.y - b[3]});
        if (std::hypot(bx, by) >= dmin) continue;
        if (closest_point(*lcurve[q], Vec3{m.x, m.y, lcurve[q]->control_points()[0].z}).distance <
            dmin)
          keep = false;
      }
      if (!keep || winding(bez_y, loops, m, eps) <= 0) continue;
      auto fc = cu[k] > c.domain_start() || cu[k + 1] < c.domain_end()
                    ? c.sub_curve(cu[k], cu[k + 1])
                    : c;
      Fragment f{fc,
                 raw[i].loop,
                 fc.evaluate(fc.domain_start()),
                 fc.evaluate(fc.domain_end()),
                 start_tangent(fc),
                 end_tangent(fc),
                 chord_length(fc),
                 false};
      if (f.length <= ctol) continue;
      frags.push_back(std::move(f));
    }
  }
  // chain into closed loops, longest fragments first
  std::vector<std::size_t> order(frags.size());
  for (std::size_t i = 0; i < order.size(); ++i) order[i] = i;
  std::sort(order.begin(), order.end(),
            [&](std::size_t x, std::size_t y) { return frags[x].length > frags[y].length; });
  std::vector<HatchContour> out;
  for (const std::size_t first : order) {
    if (frags[first].used) continue;
    std::vector<std::size_t> chain{first};
    frags[first].used = true;
    bool closed = false;
    for (;;) {
      const auto& cur = frags[chain.back()];
      if (distance(cur.b, frags[first].a) <= ctol) {
        closed = true;
        break;
      }
      std::size_t next = frags.size();
      double best = -10.0;
      for (std::size_t j = 0; j < frags.size(); ++j) {
        if (frags[j].used || distance(frags[j].a, cur.b) > ctol) continue;
        const double turn = std::atan2(cross2(cur.tb, frags[j].ta), dot2(cur.tb, frags[j].ta));
        if (turn > best) {
          best = turn;
          next = j;
        }
      }
      if (next == frags.size()) break;
      frags[next].used = true;
      chain.push_back(next);
    }
    double len = 0.0;
    for (const auto i : chain) len += frags[i].length;
    if (!closed) {
      if (len < d) continue;  // a leftover overlap at a nearly smooth corner
      throw GeometryError("hatch: contour at distance " + fmt(d) + " does not close near (" +
                          fmt(frags[chain.back()].b.x) + ", " + fmt(frags[chain.back()].b.y) + ")");
    }
    HatchContour hc;
    hc.distance = d;
    for (const auto i : chain) {
      hc.curves.push_back(frags[i].curve);
      hc.loops.push_back(frags[i].loop);
    }
    if (std::abs(loop_area(hc.curves)) <= ctol * len) continue;  // degenerate sliver
    std::sort(hc.loops.begin(), hc.loops.end());
    hc.loops.erase(std::unique(hc.loops.begin(), hc.loops.end()), hc.loops.end());
    hc.loop = hc.loops.front();
    // start the chain at the beginning of a piece of the lowest source loop
    out.push_back(std::move(hc));
  }
  std::sort(out.begin(), out.end(),
            [](const HatchContour& x, const HatchContour& y) { return x.loops < y.loops; });
  return out;
}

}  // namespace

// ------------------------------------------------------------------ hatch

Hatching hatch(const std::vector<std::vector<NurbsCurve>>& loops, const Vec3& direction,
               double spacing, const HatchOptions& options) {
  if (!(spacing > 0.0)) throw GeometryError("hatch: spacing must be > 0");
  const Vec3 d = unit_xy(direction);
  if (d.x == 0.0 && d.y == 0.0) throw GeometryError("hatch: direction has no xy component");
  const Vec3 n{-d.y, d.x, 0.0};

  Hatching out;
  if (loops.empty()) return out;

  // model size and tolerances
  double xmin = std::numeric_limits<double>::max(), xmax = -xmin, ymin = xmin, ymax = -xmin;
  for (const auto& l : loops)
    for (const auto& c : l)
      for (const auto& p : c.control_points()) {
        xmin = std::min(xmin, p.x);
        xmax = std::max(xmax, p.x);
        ymin = std::min(ymin, p.y);
        ymax = std::max(ymax, p.y);
      }
  const double size = std::max(std::hypot(xmax - xmin, ymax - ymin), 1e-300);
  const double eps = options.tolerance > 0 ? options.tolerance : 1e-9 * size;
  const double fit_tol = options.contour_tolerance > 0 ? options.contour_tolerance : 1e-7 * size;

  std::vector<double> area(loops.size());
  for (std::size_t i = 0; i < loops.size(); ++i) area[i] = loop_area(loops[i]);

  // Bezier pieces with control heights along n
  std::vector<std::vector<Bez>> bez(loops.size());
  double hmin = std::numeric_limits<double>::max(), hmax = -hmin;
  for (std::size_t li = 0; li < loops.size(); ++li) {
    if (loops[li].empty()) throw GeometryError("hatch: loop " + std::to_string(li) + " is empty");
    for (std::size_t ci = 0; ci < loops[li].size(); ++ci)
      for (const auto& b : loops[li][ci].decompose_bezier()) {
        Bez z;
        z.curve = ci;
        z.u0 = b.domain_start();
        z.u1 = b.domain_end();
        if (!(z.u1 > z.u0)) continue;
        z.w = b.weights();
        z.h.reserve(b.num_control_points());
        for (const auto& p : b.control_points()) z.h.push_back(dot2(n, p));
        z.hmin = *std::min_element(z.h.begin(), z.h.end());
        z.hmax = *std::max_element(z.h.begin(), z.h.end());
        hmin = std::min(hmin, z.hmin);
        hmax = std::max(hmax, z.hmax);
        bez[li].push_back(std::move(z));
      }
  }

  // hatch lines
  const double k0 = std::ceil((hmin - options.phase) / spacing),
               k1 = std::floor((hmax - options.phase) / spacing);
  if (k1 - k0 > 1e7) throw GeometryError("hatch: more than 10^7 hatch lines (spacing too small)");
  std::vector<Run> runs;
  std::vector<Crossing> xs;
  std::vector<int> lw(loops.size());
  for (double kk = k0; kk <= k1; kk += 1.0) {
    const auto k = static_cast<std::int64_t>(kk);
    const double c = options.phase + static_cast<double>(k) * spacing;
    xs.clear();
    for (std::size_t li = 0; li < loops.size(); ++li)
      loop_crossings(bez[li], loops[li], li, c, eps, d, runs, xs);
    if (xs.empty()) continue;
    // at equal positions, leave before entering (touching regions give separate segments)
    std::sort(xs.begin(), xs.end(), [](const Crossing& a, const Crossing& b) {
      return a.x < b.x || (a.x == b.x && a.wind < b.wind);
    });
    std::fill(lw.begin(), lw.end(), 0);
    int total = 0;
    HatchSegment cur;
    double x_start = 0.0;
    const std::size_t first_new = out.segments.size();
    for (const auto& x : xs) {
      const int before = total;
      total += x.wind;
      lw[x.loop] += x.wind;
      if (before <= 0 && total > 0) {
        x_start = x.x;
        cur = HatchSegment{};
        cur.line = k;
        cur.start_loop = x.loop;
        // region: the smallest outer loop around the segment
        double best = std::numeric_limits<double>::max();
        cur.region = x.loop;
        for (std::size_t li = 0; li < loops.size(); ++li)
          if (lw[li] != 0 && area[li] > 0 && area[li] < best) {
            best = area[li];
            cur.region = li;
          }
      } else if (before > 0 && total <= 0) {
        // merge with the previous segment of this line if they touch in the same region
        if (out.segments.size() > first_new) {
          auto& prev = out.segments.back();
          const double x_prev = dot2(d, prev.end);
          if (x_start - x_prev <= eps && prev.region == cur.region) {
            prev.end = x.x * d + c * n;
            prev.end_loop = x.loop;
            continue;
          }
        }
        cur.start = x_start * d + c * n;
        cur.end = x.x * d + c * n;
        cur.end_loop = x.loop;
        out.segments.push_back(cur);
      }
    }
    // drop zero-length and short segments
    const double lmin = std::max(options.min_length, eps);
    out.segments.erase(std::remove_if(out.segments.begin() + static_cast<std::ptrdiff_t>(first_new),
                                      out.segments.end(),
                                      [&](const HatchSegment& s) { return s.length() <= lmin; }),
                       out.segments.end());
  }

  // contours
  for (const double dist : options.contours)
    if (!(dist > 0.0)) throw GeometryError("hatch: contour distances must be > 0");
  if (!options.contours.empty()) {
    double join_gap = 0.0;  // largest gap between consecutive curves of a loop
    for (const auto& l : loops)
      for (std::size_t i = 0; i < l.size(); ++i)
        join_gap = std::max(
            join_gap,
            distance(l[i].evaluate(l[i].domain_end()),
                     l[(i + 1) % l.size()].evaluate(l[(i + 1) % l.size()].domain_start())));
    const auto bez_y = make_bez(loops, {0, 1, 0});
    for (std::size_t ci = 0; ci < options.contours.size(); ++ci)
      for (auto& hc : contours_at(loops, bez_y, options.contours[ci], fit_tol, eps, join_gap)) {
        hc.index = ci;
        out.contours.push_back(std::move(hc));
      }
    std::stable_sort(out.contours.begin(), out.contours.end(),
                     [](const HatchContour& x, const HatchContour& y) { return x.loop < y.loop; });
  }
  return out;
}

Hatching hatch(const std::vector<SectionLoop>& loops, const Vec3& direction, double spacing,
               const HatchOptions& options) {
  std::vector<std::vector<NurbsCurve>> curves(loops.size());
  for (std::size_t i = 0; i < loops.size(); ++i)
    for (const auto& sc : loops[i].curves) curves[i].push_back(sc.curve);
  return hatch(curves, direction, spacing, options);
}

}  // namespace nf
