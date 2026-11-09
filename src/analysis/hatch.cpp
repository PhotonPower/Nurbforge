#include "nurbforge/analysis/hatch.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <optional>
#include <string>

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

/// Inward offset of one loop by d (see the header).
std::vector<NurbsCurve> contour(const std::vector<NurbsCurve>& loop, std::size_t li, double d,
                                double fit_tol, double eps) {
  std::vector<NurbsCurve> pieces;
  for (const auto& c : loop)
    for (auto& p : split_at_c0_knots(c)) pieces.push_back(std::move(p));
  const std::size_t n = pieces.size();
  const auto fail = [&](const std::string& why) {
    throw GeometryError("hatch: contour at distance " + fmt(d) + " of loop " + std::to_string(li) +
                        ": " + why);
  };
  std::vector<NurbsCurve> off;
  off.reserve(n);
  for (std::size_t i = 0; i < n; ++i) {
    try {
      off.push_back(offset_curve(pieces[i], d, {0, 0, 1}, fit_tol));
    } catch (const GeometryError& e) {
      fail(std::string("offset of piece ") + std::to_string(i) + " failed (" + e.what() + ")");
    }
  }
  std::vector<double> lo(n), hi(n);
  for (std::size_t i = 0; i < n; ++i) {
    lo[i] = off[i].domain_start();
    hi[i] = off[i].domain_end();
  }
  // corner between piece i and piece i + 1: 0 smooth, arc for reflex corners
  std::vector<std::optional<NurbsCurve>> arcs(n);
  for (std::size_t i = 0; i < n; ++i) {
    const std::size_t j = (i + 1) % n;
    const Vec3 ta = end_tangent(pieces[i]), tb = start_tangent(pieces[j]);
    const double cr = cross2(ta, tb), dt = dot2(ta, tb);
    if (std::abs(cr) <= 1e-9 && dt > 0) continue;  // tangent continuous
    if (std::abs(cr) <= 1e-9) fail("the loop reverses its direction at a corner");
    const Vec3 v = pieces[i].evaluate(pieces[i].domain_end());
    if (cr > 0) {
      // convex corner (left turn): the offsets overlap, trim both at their intersection
      const auto xs = intersect(off[i], off[j], Tolerance{eps, 1e-12, 1e-9});
      const Vec3 corner = v + (d / std::cos(0.5 * std::atan2(cr, dt))) *
                                  unit_xy(Vec3{-ta.y, ta.x, 0} + Vec3{-tb.y, tb.x, 0});
      double best = std::numeric_limits<double>::max();
      double s = 0, t = 0;
      for (const auto& x : xs) {
        const double dist = distance(Vec3{x.point.x, x.point.y, 0}, Vec3{corner.x, corner.y, 0});
        if (dist < best) {
          best = dist;
          s = x.s;
          t = x.t;
        }
      }
      if (xs.empty()) fail("the offsets of the pieces at a corner do not meet (a piece vanishes)");
      hi[i] = std::min(hi[i], s);
      lo[j] = std::max(lo[j], t);
    } else {
      // reflex corner (right turn): arc about the corner from left(ta) clockwise to left(tb)
      const Vec3 x_axis{-ta.y, ta.x, 0}, y_axis{ta.x, ta.y, 0};
      const double phi = std::atan2(-cr, dt);
      arcs[i] = make_arc(Vec3{v.x, v.y, off[i].evaluate(hi[i]).z}, x_axis, y_axis, d, 0.0, phi);
    }
  }
  std::vector<NurbsCurve> out;
  for (std::size_t i = 0; i < n; ++i) {
    const double span = off[i].domain_end() - off[i].domain_start();
    if (hi[i] - lo[i] <= 1e-12 * span) fail("piece " + std::to_string(i) + " vanishes at a corner");
    out.push_back(lo[i] > off[i].domain_start() || hi[i] < off[i].domain_end()
                      ? off[i].sub_curve(lo[i], hi[i])
                      : off[i]);
    if (arcs[i]) out.push_back(*arcs[i]);
  }
  return out;
}

/// True if two pieces of contours intersect away from a shared end point.
bool pieces_cross(const NurbsCurve& a, const NurbsCurve& b, bool adjacent_ab, bool adjacent_ba,
                  double eps) {
  const auto xs = intersect(a, b, Tolerance{eps, 1e-12, 1e-9});
  const Vec3 ae = a.evaluate(a.domain_end()), as = a.evaluate(a.domain_start());
  const Vec3 bs = b.evaluate(b.domain_start()), be = b.evaluate(b.domain_end());
  for (const auto& x : xs) {
    if (adjacent_ab && distance(x.point, ae) <= 100 * eps && distance(x.point, bs) <= 100 * eps)
      continue;
    if (adjacent_ba && distance(x.point, as) <= 100 * eps && distance(x.point, be) <= 100 * eps)
      continue;
    return true;
  }
  return false;
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
  for (std::size_t li = 0; li < loops.size(); ++li)
    for (std::size_t ci = 0; ci < options.contours.size(); ++ci)
      out.contours.push_back({contour(loops[li], li, options.contours[ci], fit_tol, eps), li, ci,
                              options.contours[ci]});

  // the contours of one distance must be simple and disjoint
  for (std::size_t ci = 0; ci < options.contours.size(); ++ci) {
    std::vector<const HatchContour*> cs;
    for (const auto& hc : out.contours)
      if (hc.index == ci) cs.push_back(&hc);
    for (std::size_t a = 0; a < cs.size(); ++a)
      for (std::size_t b = a; b < cs.size(); ++b) {
        const auto& A = cs[a]->curves;
        const auto& B = cs[b]->curves;
        for (std::size_t i = 0; i < A.size(); ++i)
          for (std::size_t j = (a == b ? i + 1 : 0); j < B.size(); ++j) {
            const bool ab = a == b && j == i + 1;
            const bool ba = a == b && i == 0 && j + 1 == A.size();
            if (pieces_cross(A[i], B[j], ab, ba, fit_tol)) {
              throw GeometryError(
                  "hatch: contour at distance " + fmt(options.contours[ci]) + " of loop " +
                  std::to_string(cs[a]->loop) +
                  (a == b ? " intersects itself"
                          : " intersects the contour of loop " + std::to_string(cs[b]->loop)) +
                  " (the distance exceeds the local half-width of the region)");
            }
          }
      }
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
