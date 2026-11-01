#include "nurbforge/topology/trimmed_surface.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>
#include <optional>
#include <string>

#include "nurbforge/analysis/properties.hpp"
#include "nurbforge/analysis/surface_intersect.hpp"
#include "nurbforge/core/error.hpp"
#include "nurbforge/fitting/fit.hpp"
#include "nurbforge/geometry/primitives.hpp"
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
  poly_.push_back(outer_.polyline());
  for (const auto& h : holes_) poly_.push_back(h.polyline());
}

bool TrimmedSurface::contains(double u, double v) const {
  if (!inside_polygon(poly_[0], u, v)) return false;
  for (std::size_t i = 1; i < poly_.size(); ++i)
    if (inside_polygon(poly_[i], u, v)) return false;
  return true;
}

TrimmedSurface TrimmedSurface::with_hole(TrimLoop hole) const {
  std::vector<TrimLoop> h = holes_;
  h.push_back(std::move(hole));
  return {surface_, outer_, std::move(h)};
}

double TrimmedSurface::area(double rel_tol) const {
  // Green's theorem in the parameter plane: with F(u, v) = int_{u0}^{u} |S_u x S_v|(s, v) ds,
  // area = sum over the loops of int F(u(t), v(t)) v'(t) dt (holes are clockwise: negative).
  const double U0 = surface_.domain_u()[0];
  const auto uknots = unique_knots(surface_.knots_u()).first;
  const auto dA = [&](double u, double v) {
    const auto d = surface_.first_derivatives(u, v);
    return norm(cross(d[1], d[2]));
  };
  const auto F = [&](double u, double v) {
    return integrate([&](double s) { return dA(s, v); }, U0, u, 0.1 * rel_tol, uknots);
  };
  // absolute error target relative to the untrimmed area: loop integrands may be pure rounding
  // noise (curves along u, where v' ~ 0), so a relative target alone would never be met
  const double abs_tol = rel_tol * nf::area(surface_, 1e-6);
  double total = 0.0;
  const auto loop_integral = [&](const TrimLoop& l) {
    double sum = 0.0;
    for (const auto& c : l.curves) {
      std::vector<double> br = unique_knots(c.knots()).first;
      sum += integrate(
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
    // piecewise per span of the trim curve: cubic interpolation of S(pc(t)) at the original
    // parameter, doubled until the midpoint error is below the tolerance
    std::optional<NurbsCurve> result;
    const auto br = unique_knots(pc.knots()).first;
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
