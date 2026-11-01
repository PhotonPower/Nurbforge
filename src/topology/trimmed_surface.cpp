#include "nurbforge/topology/trimmed_surface.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <optional>
#include <string>

#include "nurbforge/analysis/closest_point.hpp"
#include "nurbforge/analysis/properties.hpp"
#include "nurbforge/analysis/surface_intersect.hpp"
#include "nurbforge/core/error.hpp"
#include "nurbforge/fitting/fit.hpp"
#include "nurbforge/geometry/primitives.hpp"
#include "nurbforge/geometry/transform.hpp"
#include "nurbforge/nurbs/join.hpp"

namespace nf {

namespace {

using P2 = std::array<double, 2>;

NurbsCurve line2(const P2& a, const P2& b) { return make_line({a[0], a[1], 0.0}, {b[0], b[1], 0.0}); }

P2 start2(const NurbsCurve& c) {
  const Vec3 p = c.evaluate(c.domain_start());
  return {p.x, p.y};
}
P2 end2(const NurbsCurve& c) {
  const Vec3 p = c.evaluate(c.domain_end());
  return {p.x, p.y};
}
double dist2(const P2& a, const P2& b) { return std::hypot(a[0] - b[0], a[1] - b[1]); }

// Even-odd point in polygon test.
bool inside_polygon(const std::vector<P2>& poly, double u, double v) {
  bool in = false;
  for (std::size_t i = 0, j = poly.size() - 1; i < poly.size(); j = i++) {
    const auto& a = poly[i];
    const auto& b = poly[j];
    if ((a[1] > v) != (b[1] > v) && u < (b[0] - a[0]) * (v - a[1]) / (b[1] - a[1]) + a[0]) in = !in;
  }
  return in;
}

double shoelace(const std::vector<P2>& poly) {
  double s = 0.0;
  for (std::size_t i = 0, j = poly.size() - 1; i < poly.size(); j = i++)
    s += (poly[j][0] - poly[i][0]) * (poly[j][1] + poly[i][1]);
  return 0.5 * s;
}

}  // namespace

// ------------------------------------------------------------------ TrimLoop

std::vector<P2> TrimLoop::polyline(std::size_t samples_per_span) const {
  std::vector<P2> pts;
  for (const auto& c : curves) {
    const std::size_t spans = std::max<std::size_t>(1, unique_knots(c.knots()).first.size() - 1);
    const std::size_t n = spans * std::max<std::size_t>(samples_per_span, 2);
    const double a = c.domain_start(), b = c.domain_end();
    for (std::size_t i = 0; i < n; ++i) {  // last point omitted: it starts the next curve
      const Vec3 p = c.evaluate(a + (b - a) * static_cast<double>(i) / static_cast<double>(n));
      pts.push_back({p.x, p.y});
    }
  }
  return pts;
}

double TrimLoop::signed_area() const { return shoelace(polyline()); }

TrimLoop TrimLoop::reversed() const {
  TrimLoop r;
  for (auto it = curves.rbegin(); it != curves.rend(); ++it) r.curves.push_back(it->reversed());
  return r;
}

TrimLoop rectangle_loop(double u0, double u1, double v0, double v1) {
  if (!(u1 > u0) || !(v1 > v0)) throw GeometryError("rectangle_loop: empty rectangle");
  return {{line2({u0, v0}, {u1, v0}), line2({u1, v0}, {u1, v1}), line2({u1, v1}, {u0, v1}),
           line2({u0, v1}, {u0, v0})}};
}

TrimLoop circle_loop(double cu, double cv, double r) {
  return {{make_circle({cu, cv, 0.0}, {0, 0, 1}, r)}};
}

TrimLoop polygon_loop(const std::vector<P2>& points) {
  if (points.size() < 3) throw GeometryError("polygon_loop: need at least three points");
  TrimLoop l;
  for (std::size_t i = 0; i < points.size(); ++i) l.curves.push_back(line2(points[i], points[(i + 1) % points.size()]));
  return l;
}

// ------------------------------------------------------------------ TrimmedSurface

TrimmedSurface::TrimmedSurface(NurbsSurface surface) : surface_(std::move(surface)) {
  const auto du = surface_.domain_u(), dv = surface_.domain_v();
  outer_ = rectangle_loop(du[0], du[1], dv[0], dv[1]);
  build_polylines();
}

TrimmedSurface::TrimmedSurface(NurbsSurface surface, TrimLoop outer, std::vector<TrimLoop> holes,
                               const Tolerance& tol)
    : surface_(std::move(surface)), outer_(std::move(outer)), holes_(std::move(holes)) {
  const auto du = surface_.domain_u(), dv = surface_.domain_v();
  const double gap_tol = std::max(tol.length, 1e-9 * std::max(du[1] - du[0], dv[1] - dv[0]));
  const auto check = [&](const TrimLoop& l) {
    if (l.curves.empty()) throw GeometryError("trimmed surface: empty loop");
    for (std::size_t i = 0; i < l.curves.size(); ++i) {
      const P2 e = end2(l.curves[i]), s = start2(l.curves[(i + 1) % l.curves.size()]);
      if (dist2(e, s) > gap_tol) throw GeometryError("trimmed surface: loop is not closed");
    }
    const double eu = 1e-9 * (du[1] - du[0]), ev = 1e-9 * (dv[1] - dv[0]);
    for (const auto& p : l.polyline(4))
      if (p[0] < du[0] - eu || p[0] > du[1] + eu || p[1] < dv[0] - ev || p[1] > dv[1] + ev)
        throw GeometryError("trimmed surface: loop leaves the parameter domain");
  };
  check(outer_);
  for (const auto& h : holes_) check(h);
  if (outer_.signed_area() < 0.0) outer_ = outer_.reversed();
  for (auto& h : holes_)
    if (h.signed_area() > 0.0) h = h.reversed();
  build_polylines();
}

void TrimmedSurface::build_polylines() {
  poly_.clear();
  band_.clear();
  const auto du = surface_.domain_u(), dv = surface_.domain_v();
  const double tiny = 1e-12 * std::max(du[1] - du[0], dv[1] - dv[0]);
  const auto add = [&](const TrimLoop& l) {
    poly_.push_back(l.polyline());
    // deviation of every curve from its polyline chords, at the chord midpoints (same sampling
    // as polyline()): inside this band the polygon test is not reliable
    double dev = 0.0;
    for (const auto& c : l.curves) {
      const std::size_t spans = std::max<std::size_t>(1, unique_knots(c.knots()).first.size() - 1);
      const std::size_t n = spans * 16;
      const double a = c.domain_start(), b = c.domain_end();
      for (std::size_t i = 0; i < n; ++i) {
        const auto at = [&](double f) {
          const Vec3 p = c.evaluate(a + (b - a) * f);
          return P2{p.x, p.y};
        };
        const P2 p0 = at(static_cast<double>(i) / static_cast<double>(n));
        const P2 p1 = at(static_cast<double>(i + 1) / static_cast<double>(n));
        const P2 pm = at((static_cast<double>(i) + 0.5) / static_cast<double>(n));
        const double ex = p1[0] - p0[0], ey = p1[1] - p0[1], len = std::hypot(ex, ey);
        if (len > 0.0) dev = std::max(dev, std::abs((pm[0] - p0[0]) * ey - (pm[1] - p0[1]) * ex) / len);
      }
    }
    band_.push_back(2.0 * dev + tiny);
  };
  add(outer_);
  for (const auto& h : holes_) add(h);
}

// Exact side test near a loop: the closest point on the loop's curves and the curve direction
// there (loops are oriented with the face on their left: outer counter-clockwise, holes
// clockwise). At a corner (foot at a curve end) the point is on the face side of a convex corner
// if it is left of both adjacent curves, of a reflex corner if it is left of either.
bool TrimmedSurface::left_of_loop(const TrimLoop& loop, double u, double v) const {
  const Vec3 P{u, v, 0.0};
  std::size_t best = 0;
  CurvePoint foot;
  foot.distance = std::numeric_limits<double>::max();
  for (std::size_t k = 0; k < loop.curves.size(); ++k) {
    const CurvePoint f = closest_point(loop.curves[k], P);
    if (f.distance < foot.distance) {
      foot = f;
      best = k;
    }
  }
  const NurbsCurve& c = loop.curves[best];
  const auto left = [&](const NurbsCurve& cu, double t) {
    const auto d = cu.derivatives(t, 1);
    return d[1].x * (P.y - d[0].y) - d[1].y * (P.x - d[0].x) > 0.0;
  };
  const double span = c.domain_end() - c.domain_start(), eps = 1e-9 * span;
  const std::size_t m = loop.curves.size();
  const NurbsCurve* in = nullptr;
  const NurbsCurve* out = nullptr;
  double t_in = 0.0, t_out = 0.0;
  if (foot.u >= c.domain_end() - eps) {
    in = &c;
    t_in = c.domain_end();
    out = &loop.curves[(best + 1) % m];
    t_out = out->domain_start();
  } else if (foot.u <= c.domain_start() + eps) {
    in = &loop.curves[(best + m - 1) % m];
    t_in = in->domain_end();
    out = &c;
    t_out = c.domain_start();
  }
  if (!in) return left(c, foot.u);
  const Vec3 ti = in->derivatives(t_in, 1)[1], to = out->derivatives(t_out, 1)[1];
  const bool convex = ti.x * to.y - ti.y * to.x >= 0.0;
  const bool l_in = left(*in, t_in), l_out = left(*out, t_out);
  return convex ? (l_in && l_out) : (l_in || l_out);
}

bool TrimmedSurface::is_untrimmed() const {
  if (!holes_.empty()) return false;
  const auto du = surface_.domain_u(), dv = surface_.domain_v();
  const double eu = 1e-9 * (du[1] - du[0]), ev = 1e-9 * (dv[1] - dv[0]);
  for (const auto& p : poly_[0]) {
    const bool on = std::abs(p[0] - du[0]) <= eu || std::abs(p[0] - du[1]) <= eu ||
                    std::abs(p[1] - dv[0]) <= ev || std::abs(p[1] - dv[1]) <= ev;
    if (!on) return false;
  }
  return std::abs(std::abs(shoelace(poly_[0])) - (du[1] - du[0]) * (dv[1] - dv[0])) <=
         1e-9 * (du[1] - du[0]) * (dv[1] - dv[0]);
}

bool TrimmedSurface::contains(double u, double v) const {
  // polygon test, refined by the exact side test within the chord-deviation band of a loop
  // (a 64-gon of a trim circle deviates by 6e-4 – points of intersection curves near the
  // boundary would otherwise be misclassified and cuts would end short of the loop)
  const auto face_side = [&](std::size_t k) {
    const TrimLoop& loop = k == 0 ? outer_ : holes_[k - 1];
    const auto& poly = poly_[k];
    double d = std::numeric_limits<double>::max();
    for (std::size_t i = 0, j = poly.size() - 1; i < poly.size(); j = i++) {
      const double ax = poly[j][0], ay = poly[j][1], ex = poly[i][0] - ax, ey = poly[i][1] - ay;
      const double l2 = ex * ex + ey * ey;
      const double t = l2 > 0.0 ? std::clamp(((u - ax) * ex + (v - ay) * ey) / l2, 0.0, 1.0) : 0.0;
      d = std::min(d, std::hypot(u - ax - t * ex, v - ay - t * ey));
    }
    if (d <= band_[k]) return left_of_loop(loop, u, v);
    const bool in = inside_polygon(poly, u, v);
    return k == 0 ? in : !in;
  };
  for (std::size_t k = 0; k < poly_.size(); ++k)
    if (!face_side(k)) return false;
  return true;
}

TrimmedSurface TrimmedSurface::flipped() const {
  const auto du = surface_.domain_u();
  const double a = du[0], b = du[1];
  // reversed u direction: rows in opposite order, mirrored knots
  const KnotVector& K = surface_.knots_u();
  KnotVector U(K.size());
  for (std::size_t i = 0; i < K.size(); ++i) U[i] = a + b - K[K.size() - 1 - i];
  std::vector<Vec3> P;
  std::vector<double> W;
  for (std::size_t i = surface_.num_u(); i-- > 0;)
    for (std::size_t j = 0; j < surface_.num_v(); ++j) {
      P.push_back(surface_.control_point(i, j));
      W.push_back(surface_.weight(i, j));
    }
  NurbsSurface s{surface_.degree_u(), surface_.degree_v(), std::move(U), surface_.knots_v(),
                 surface_.num_u(), surface_.num_v(), std::move(P), std::move(W)};
  // mirror the loops: x -> a + b - x (this reverses their orientation; the constructor fixes it)
  const Transform mirror({-1.0, 0.0, 0.0, a + b, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0});
  const auto map_loop = [&](const TrimLoop& l) {
    TrimLoop m;
    for (const auto& c : l.curves) m.curves.push_back(c.transformed(mirror));
    return m;
  };
  std::vector<TrimLoop> holes;
  for (const auto& h : holes_) holes.push_back(map_loop(h));
  return {std::move(s), map_loop(outer_), std::move(holes)};
}

TrimmedSurface TrimmedSurface::with_hole(TrimLoop hole) const {
  std::vector<TrimLoop> h = holes_;
  h.push_back(std::move(hole));
  return {surface_, outer_, std::move(h)};
}

double TrimmedSurface::area(double rel_tol) const {
  const auto dA = [&](double u, double v) {
    const auto d = surface_.first_derivatives(u, v);
    return norm(cross(d[1], d[2]));
  };
  // absolute error target relative to the untrimmed area: loop integrands may be pure rounding
  // noise (curves along u, where v' ~ 0), so a relative target alone would never be met
  return integrate(dA, rel_tol, rel_tol * nf::area(surface_, 1e-6));
}

double TrimmedSurface::integrate(const std::function<double(double, double)>& f, double rel_tol,
                                 double abs_tol) const {
  // Green's theorem in the parameter plane: with F(u, v) = int_{u0}^{u} f(s, v) ds the integral
  // is the sum over the loops of int F(u(t), v(t)) v'(t) dt (holes are clockwise: negative).
  const double U0 = surface_.domain_u()[0];
  const auto uknots = unique_knots(surface_.knots_u()).first;
  const auto F = [&](double u, double v) {
    return nf::integrate([&](double s) { return f(s, v); }, U0, u, 0.1 * rel_tol, uknots,
                         0.01 * abs_tol);
  };
  double total = 0.0;
  const auto loop_integral = [&](const TrimLoop& l) {
    double sum = 0.0;
    for (const auto& c : l.curves) {
      // break only where the curve is less smooth than C2 (multiplicity >= p - 1): interpolated
      // intersection curves have hundreds of simple (C2) knots, which need no breaks
      std::vector<double> br;
      {
        const auto [vals, mults] = unique_knots(c.knots());
        const int need = std::max(1, static_cast<int>(c.degree()) - 1);
        for (std::size_t q = 1; q + 1 < vals.size(); ++q)
          if (mults[q] >= need) br.push_back(vals[q]);
      }
      sum += nf::integrate(
          [&](double t) {
            const auto d = c.derivatives(t, 1);
            return d[1].y == 0.0 ? 0.0 : F(d[0].x, d[0].y) * d[1].y;
          },
          c.domain_start(), c.domain_end(), rel_tol, br, abs_tol / static_cast<double>(4 * l.curves.size()));
    }
    return sum;
  };
  total += loop_integral(outer_);
  for (const auto& h : holes_) total += loop_integral(h);
  return total;
}

std::vector<std::vector<NurbsCurve>> TrimmedSurface::boundary_curves(double tolerance) const {
  std::vector<std::vector<NurbsCurve>> out;
  const auto map_curve = [&](const NurbsCurve& pc) {
    // piecewise per span of the trim curve and between crossings of the surface's knot lines
    // (where S is only C1, e.g. rational circles – a cubic across them converges like h^2):
    // cubic interpolation of S(pc(t)) at the original parameter, doubled until the midpoint
    // error is below the tolerance
    std::optional<NurbsCurve> result;
    std::vector<double> br = unique_knots(pc.knots()).first;
    {
      const auto ku = unique_knots(surface_.knots_u()).first, kv = unique_knots(surface_.knots_v()).first;
      std::vector<double> extra;
      // A crossing needs the curve to be clearly off the knot line on both sides (|g| > eps):
      // a trim curve running along a knot line (e.g. an intersection with a plane through a
      // sphere's meridian) must not produce a crossing at every rounding-level sign change.
      const auto crossings = [&](int coord, const std::vector<double>& knots, double eps) {
        for (std::size_t s = 0; s + 1 < br.size(); ++s)
          for (double kappa : knots) {
            const auto g = [&](double t) {
              const Vec3 w = pc.evaluate(t);
              const double x = (coord == 0 ? w.x : w.y) - kappa;
              return std::abs(x) <= eps ? 0.0 : x;
            };
            const int n = 32;
            double t0 = br[s], g0 = g(t0);
            for (int i = 1; i <= n; ++i) {
              const double t1 = br[s] + (br[s + 1] - br[s]) * i / n, g1 = g(t1);
              if (g1 == 0.0 && g0 != 0.0 && i < n) {  // knot value hit exactly at a sample
                extra.push_back(t1);
              } else if (g0 * g1 < 0.0) {  // strict sign change: bisection
                double lo = t0, hi = t1;
                for (int it = 0; it < 60; ++it) {
                  const double m = 0.5 * (lo + hi);
                  (g(m) * g0 <= 0.0 ? hi : lo) = m;
                }
                extra.push_back(0.5 * (lo + hi));
              }
              t0 = t1;
              g0 = g1;
            }
          }
      };
      crossings(0, std::vector<double>(ku.begin() + 1, ku.end() - 1), 1e-9 * (ku.back() - ku.front()));
      crossings(1, std::vector<double>(kv.begin() + 1, kv.end() - 1), 1e-9 * (kv.back() - kv.front()));
      br.insert(br.end(), extra.begin(), extra.end());
      std::sort(br.begin(), br.end());
      br.erase(std::unique(br.begin(), br.end(), [](double x, double y) { return y - x <= 1e-12; }), br.end());
    }
    for (std::size_t s = 0; s + 1 < br.size(); ++s) {
      const double a = br[s], b = br[s + 1];
      for (std::size_t N = 4;; N = 2 * N - 1) {
        std::vector<double> t(N);
        std::vector<Vec3> P(N);
        for (std::size_t j = 0; j < N; ++j) {
          t[j] = j + 1 == N ? b : a + (b - a) * static_cast<double>(j) / static_cast<double>(N - 1);
          const Vec3 w = pc.evaluate(t[j]);
          P[j] = surface_.evaluate(w.x, w.y);
        }
        const NurbsCurve piece = interpolate_curve(P, 3, t);
        double err = 0.0;
        for (std::size_t j = 0; j + 1 < N; ++j) {
          const double tm = 0.5 * (t[j] + t[j + 1]);
          const Vec3 w = pc.evaluate(tm);
          err = std::max(err, distance(piece.evaluate(tm), surface_.evaluate(w.x, w.y)));
        }
        if (err <= tolerance || N > 8193) {
          result = result ? join(*result, piece, 10 * tolerance) : piece;
          break;
        }
      }
    }
    return *result;
  };
  for (const TrimLoop* l : [&] {
         std::vector<const TrimLoop*> v{&outer_};
         for (const auto& h : holes_) v.push_back(&h);
         return v;
       }()) {
    std::vector<NurbsCurve> loop;
    for (const auto& c : l->curves) loop.push_back(map_curve(c));
    out.push_back(std::move(loop));
  }
  return out;
}

// ------------------------------------------------------------------ split

std::vector<TrimmedSurface> split(const NurbsSurface& surface, const std::vector<NurbsCurve>& pcurves,
                                  const Tolerance& tol) {
  const auto du = surface.domain_u(), dv = surface.domain_v();
  const double U0 = du[0], U1 = du[1], V0 = dv[0], V1 = dv[1];
  const double eps = std::max(10 * tol.parametric, 1e-9 * std::max(U1 - U0, V1 - V0));
  const auto on_boundary = [&](const P2& p) {
    return std::abs(p[0] - U0) <= eps || std::abs(p[0] - U1) <= eps || std::abs(p[1] - V0) <= eps ||
           std::abs(p[1] - V1) <= eps;
  };
  // perimeter coordinate, counter-clockwise: bottom [0,1), right [1,2), top [2,3), left [3,4)
  const auto perim = [&](const P2& p) {
    if (std::abs(p[1] - V0) <= eps && p[0] < U1 - eps) return (p[0] - U0) / (U1 - U0);
    if (std::abs(p[0] - U1) <= eps && p[1] < V1 - eps) return 1.0 + (p[1] - V0) / (V1 - V0);
    if (std::abs(p[1] - V1) <= eps && p[0] > U0 + eps) return 2.0 + (U1 - p[0]) / (U1 - U0);
    return 3.0 + (V1 - p[1]) / (V1 - V0);
  };
  const auto at = [&](double s) -> P2 {
    s = std::fmod(s, 4.0);
    if (s < 1.0) return {U0 + s * (U1 - U0), V0};
    if (s < 2.0) return {U1, V0 + (s - 1.0) * (V1 - V0)};
    if (s < 3.0) return {U1 - (s - 2.0) * (U1 - U0), V1};
    return {U0, V1 - (s - 3.0) * (V1 - V0)};
  };
  // boundary lines from perimeter position s0 to s1 (counter-clockwise), split at corners
  const auto boundary = [&](double s0, double s1, std::vector<NurbsCurve>& out) {
    if (s1 < s0) s1 += 4.0;
    double s = s0;
    while (s1 - s > 1e-12) {
      const double next = std::min(s1, std::floor(s + 1e-12) + 1.0);
      const P2 a = at(s), b = at(next);
      if (dist2(a, b) > eps) out.push_back(line2(a, b));
      s = next;
    }
  };

  std::vector<NurbsCurve> chords;
  std::vector<TrimLoop> closed;
  for (const auto& c : pcurves) {
    const P2 a = start2(c), b = end2(c);
    if (dist2(a, b) <= eps) {
      closed.push_back({{c}});
    } else if (on_boundary(a) && on_boundary(b)) {
      chords.push_back(c);
    } else {
      throw GeometryError("split: curve ends inside the parameter domain (dangling curve)");
    }
  }

  // regions of the domain cut by the chords (non-crossing chords: #chords + 1 regions)
  std::vector<TrimLoop> regions;
  if (chords.empty()) {
    regions.push_back(rectangle_loop(U0, U1, V0, V1));
  } else {
    struct End {
      double s;
      std::size_t chord;
      bool at_start;
    };
    std::vector<End> ends;
    for (std::size_t i = 0; i < chords.size(); ++i) {
      ends.push_back({perim(start2(chords[i])), i, true});
      ends.push_back({perim(end2(chords[i])), i, false});
    }
    std::sort(ends.begin(), ends.end(), [](const End& x, const End& y) { return x.s < y.s; });
    const std::size_t N = ends.size();
    const auto partner = [&](std::size_t k) {
      for (std::size_t m = 0; m < N; ++m)
        if (m != k && ends[m].chord == ends[k].chord) return m;
      return k;
    };
    std::vector<bool> used(N, false);  // boundary arc from ends[k] to ends[k+1]
    for (std::size_t k0 = 0; k0 < N; ++k0) {
      if (used[k0]) continue;
      TrimLoop loop;
      std::size_t k = k0;
      for (std::size_t guard = 0; guard <= N; ++guard) {
        used[k] = true;
        const std::size_t next = (k + 1) % N;
        boundary(ends[k].s, ends[next].s + (next == 0 ? 4.0 : 0.0), loop.curves);
        const NurbsCurve& ch = chords[ends[next].chord];
        loop.curves.push_back(ends[next].at_start ? ch : ch.reversed());
        k = partner(next);
        if (k == k0) break;
      }
      regions.push_back(std::move(loop));
    }
  }

  // closed loops: parent = innermost closed loop containing them, else the region around them
  const std::size_t R = regions.size(), C = closed.size();
  std::vector<std::vector<P2>> rpoly(R), cpoly(C);
  std::vector<double> carea(C);
  for (std::size_t r = 0; r < R; ++r) rpoly[r] = regions[r].polyline();
  for (std::size_t i = 0; i < C; ++i) {
    if (closed[i].signed_area() < 0.0) closed[i] = closed[i].reversed();
    cpoly[i] = closed[i].polyline();
    carea[i] = std::abs(shoelace(cpoly[i]));
  }
  std::vector<std::ptrdiff_t> parent(C, -1);  // index into closed, or -(r+2) for region r
  for (std::size_t i = 0; i < C; ++i) {
    const P2 p = cpoly[i].front();
    double best = 1e300;
    for (std::size_t j = 0; j < C; ++j)
      if (j != i && carea[j] > carea[i] && carea[j] < best && inside_polygon(cpoly[j], p[0], p[1])) {
        best = carea[j];
        parent[i] = static_cast<std::ptrdiff_t>(j);
      }
    if (parent[i] == -1)
      for (std::size_t r = 0; r < R; ++r)
        if (inside_polygon(rpoly[r], p[0], p[1])) parent[i] = -static_cast<std::ptrdiff_t>(r) - 2;
  }
  std::vector<TrimmedSurface> out;
  for (std::size_t r = 0; r < R; ++r) {
    std::vector<TrimLoop> holes;
    for (std::size_t i = 0; i < C; ++i)
      if (parent[i] == -static_cast<std::ptrdiff_t>(r) - 2) holes.push_back(closed[i]);
    out.emplace_back(surface, regions[r], std::move(holes), tol);
  }
  for (std::size_t i = 0; i < C; ++i) {
    std::vector<TrimLoop> holes;
    for (std::size_t j = 0; j < C; ++j)
      if (parent[j] == static_cast<std::ptrdiff_t>(i)) holes.push_back(closed[j]);
    out.emplace_back(surface, closed[i], std::move(holes), tol);
  }
  return out;
}

std::vector<TrimmedSurface> split(const NurbsSurface& surface, const NurbsSurface& cutter,
                                  const Tolerance& tol) {
  SurfaceIntersectionOptions opt;
  opt.tol = tol;
  const auto res = intersect(surface, cutter, opt);
  std::vector<NurbsCurve> pcurves;
  for (const auto& c : res.curves) pcurves.push_back(c.pcurve1);
  return split(surface, pcurves, tol);
}

}  // namespace nf
