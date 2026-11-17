#include "nurbforge/topology/region.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <utility>

#include "loops.hpp"
#include "nurbforge/analysis/properties.hpp"
#include "nurbforge/core/error.hpp"
#include "nurbforge/geometry/primitives.hpp"
#include "nurbforge/nurbs/knot_vector.hpp"
#include "nurbforge/topology/trimmed_surface.hpp"

namespace nf {

namespace {

using P2 = std::array<double, 2>;

NurbsCurve map_xy(const NurbsCurve& c, double x0, double y0, double sx, double sy) {
  std::vector<Vec3> pts;
  for (const auto& p : c.control_points()) pts.push_back({x0 + sx * p.x, y0 + sy * p.y, 0.0});
  return {c.degree(), c.knots(), std::move(pts), c.weights()};
}

bool inside(const std::vector<P2>& poly, const P2& q) {
  bool in = false;
  for (std::size_t i = 0, j = poly.size() - 1; i < poly.size(); j = i++)
    if ((poly[i][1] > q[1]) != (poly[j][1] > q[1]) &&
        q[0] < poly[j][0] +
                   (poly[i][0] - poly[j][0]) * (q[1] - poly[j][1]) / (poly[i][1] - poly[j][1]))
      in = !in;
  return in;
}

// Loops (as TrimLoops, z = 0) with their nesting depth and immediate parent.
struct Nested {
  std::vector<TrimLoop> loops;
  std::vector<std::size_t> depth;
  std::vector<std::ptrdiff_t> parent;
};

Nested nest(std::vector<TrimLoop> loops) {
  Nested n;
  const std::size_t m = loops.size();
  std::vector<std::vector<P2>> poly(m);
  std::vector<double> area(m);
  for (std::size_t i = 0; i < m; ++i) {
    poly[i] = loops[i].polyline(16);
    area[i] = std::abs(loops[i].signed_area());
  }
  n.depth.assign(m, 0);
  n.parent.assign(m, -1);
  for (std::size_t i = 0; i < m; ++i) {
    double best = 1e300;
    for (std::size_t j = 0; j < m; ++j) {
      if (i == j || area[j] <= area[i] || !inside(poly[j], poly[i].front())) continue;
      ++n.depth[i];
      if (area[j] < best) {
        best = area[j];
        n.parent[i] = static_cast<std::ptrdiff_t>(j);
      }
    }
  }
  // orientation by depth: outer counter-clockwise, holes clockwise
  for (std::size_t i = 0; i < m; ++i) {
    const bool ccw = loops[i].signed_area() > 0.0;
    if (ccw != (n.depth[i] % 2 == 0)) loops[i] = loops[i].reversed();
  }
  n.loops = std::move(loops);
  return n;
}

TrimLoop as_loop(const std::vector<NurbsCurve>& curves) {
  TrimLoop l;
  for (const auto& c : curves) {
    std::vector<Vec3> pts = c.control_points();
    for (auto& p : pts) p.z = 0.0;
    l.curves.emplace_back(c.degree(), c.knots(), std::move(pts), c.weights());
  }
  return l;
}

// The faces of a region on the rectangle surface (parameters = mapped xy).
std::vector<TrimmedSurface> faces(const NurbsSurface& rect, std::vector<TrimLoop> loops) {
  const Nested n = nest(std::move(loops));
  std::vector<TrimmedSurface> out;
  for (std::size_t i = 0; i < n.loops.size(); ++i) {
    if (n.depth[i] % 2 != 0) continue;
    std::vector<TrimLoop> holes;
    for (std::size_t j = 0; j < n.loops.size(); ++j)
      if (n.depth[j] % 2 != 0 && n.parent[j] == static_cast<std::ptrdiff_t>(i))
        holes.push_back(n.loops[j]);
    out.emplace_back(rect, n.loops[i], std::move(holes));
  }
  return out;
}

// A point well inside the piece (parameter space): the grid point farthest from its loops.
P2 interior_point(const TrimmedSurface& f) {
  double u0 = 1e300, u1 = -1e300, v0 = 1e300, v1 = -1e300;
  for (const auto& p : f.outer().polyline(8)) {
    u0 = std::min(u0, p[0]), u1 = std::max(u1, p[0]);
    v0 = std::min(v0, p[1]), v1 = std::max(v1, p[1]);
  }
  std::vector<std::vector<P2>> polys{f.outer().polyline(16)};
  for (const auto& h : f.holes()) polys.push_back(h.polyline(16));
  P2 best{0.5 * (u0 + u1), 0.5 * (v0 + v1)};
  double best_d = -1.0;
  for (int n : {16, 64}) {
    for (int i = 1; i < n; ++i)
      for (int j = 1; j < n; ++j) {
        const double u = u0 + (u1 - u0) * (i + 0.37) / n, v = v0 + (v1 - v0) * (j + 0.61) / n;
        if (!f.contains(u, v)) continue;
        double d = 1e300;
        for (const auto& poly : polys)
          for (std::size_t a = 0, b = poly.size() - 1; a < poly.size(); b = a++) {
            const double ex = poly[a][0] - poly[b][0], ey = poly[a][1] - poly[b][1];
            const double l2 = ex * ex + ey * ey;
            const double t =
                l2 > 0 ? std::clamp(((u - poly[b][0]) * ex + (v - poly[b][1]) * ey) / l2, 0.0, 1.0)
                       : 0.0;
            d = std::min(d, std::hypot(u - poly[b][0] - t * ex, v - poly[b][1] - t * ey));
          }
        if (d > best_d) {
          best_d = d;
          best = {u, v};
        }
      }
    if (best_d > 0.0) break;
  }
  if (best_d <= 0.0)
    throw GeometryError("region_boolean: cannot find an interior point of a piece");
  return best;
}

}  // namespace

Region oriented(const Region& region) {
  std::vector<TrimLoop> loops;
  for (const auto& l : region) loops.push_back(as_loop(l));
  const Nested n = nest(std::move(loops));
  Region out;
  for (const auto& l : n.loops) out.push_back(l.curves);
  return out;
}

double region_area(const Region& region) {
  // Green's theorem, A = 1/2 sum of int (x y' - y x') dt over the exact curves
  double a = 0.0;
  for (const auto& l : oriented(region))
    for (const auto& c : l) {
      const auto f = [&](double t) {
        const auto d = c.derivatives(t, 1);
        return 0.5 * (d[0].x * d[1].y - d[0].y * d[1].x);
      };
      std::vector<double> breaks;
      for (const double k : unique_knots(c.knots()).first) breaks.push_back(k);
      a += integrate(f, c.domain_start(), c.domain_end(), 1e-12, breaks, 1e-15);
    }
  return a;
}

Region region_boolean(const Region& a, const Region& b, BooleanOp op, double tolerance) {
  // a rectangle covering both regions; its parameters are the affinely mapped xy coordinates
  double x0 = 1e300, x1 = -1e300, y0 = 1e300, y1 = -1e300;
  for (const Region* r : {&a, &b})
    for (const auto& l : *r)
      for (const auto& c : l)
        for (const auto& p : c.control_points()) {
          x0 = std::min(x0, p.x), x1 = std::max(x1, p.x);
          y0 = std::min(y0, p.y), y1 = std::max(y1, p.y);
        }
  if (!(x1 > x0) || !(y1 > y0)) throw GeometryError("region_boolean: empty or degenerate regions");
  const double pad = 0.05 * std::max(x1 - x0, y1 - y0);
  x0 -= pad, x1 += pad, y0 -= pad, y1 += pad;
  const double W = x1 - x0, H = y1 - y0;
  const NurbsSurface rect = make_rectangle({x0, y0, 0}, {W, 0, 0}, {0, H, 0});
  const auto to_uv = [&](const Region& r) {
    std::vector<TrimLoop> loops;
    for (const auto& l : r) {
      TrimLoop t;
      for (const auto& c : l) t.curves.push_back(map_xy(c, -x0 / W, -y0 / H, 1.0 / W, 1.0 / H));
      loops.push_back(std::move(t));
    }
    return loops;
  };
  const auto la = to_uv(a), lb = to_uv(b);
  std::vector<NurbsCurve> cuts_a, cuts_b;  // the boundary curves of each region as cuts
  for (const auto& l : la) cuts_a.insert(cuts_a.end(), l.curves.begin(), l.curves.end());
  for (const auto& l : lb) cuts_b.insert(cuts_b.end(), l.curves.begin(), l.curves.end());
  const auto fa = faces(rect, la), fb = faces(rect, lb);
  const Tolerance tol{tolerance / std::max(W, H), 1e-12, 1e-9};
  const auto in = [](const std::vector<TrimmedSurface>& fs, const P2& q) {
    for (const auto& f : fs)
      if (f.contains(q[0], q[1])) return true;
    return false;
  };
  std::vector<NurbsCurve> kept;
  const auto take =
      [&](const std::vector<TrimmedSurface>& own, const std::vector<NurbsCurve>& other_cuts,
          const std::vector<TrimmedSurface>& other, bool keep_inside, bool keep_outside) {
        for (const auto& f : own)
          for (const auto& piece : split(f, other_cuts, tol)) {
            const bool inside_other = in(other, interior_point(piece));
            if ((inside_other && keep_inside) || (!inside_other && keep_outside)) {
              for (const auto& c : piece.outer().curves) kept.push_back(c);
              for (const auto& h : piece.holes())
                for (const auto& c : h.curves) kept.push_back(c);
            }
          }
      };
  switch (op) {
    case BooleanOp::Union:
      take(fa, cuts_b, fb, true, true);
      take(fb, cuts_a, fa, false, true);
      break;
    case BooleanOp::Intersection:
      take(fa, cuts_b, fb, true, false);
      break;
    case BooleanOp::Difference:
      take(fa, cuts_b, fb, false, true);
      break;
  }
  if (kept.empty()) return {};
  const auto loops = detail::union_loops(std::move(kept), tolerance / std::max(W, H));
  if (loops.empty()) throw GeometryError("region_boolean: the result does not form closed loops");
  Region out;
  for (const auto& l : loops) {
    std::vector<NurbsCurve> curves;
    for (const auto& c : l.curves) curves.push_back(map_xy(c, x0, y0, W, H));
    out.push_back(std::move(curves));
  }
  return oriented(out);
}

}  // namespace nf
