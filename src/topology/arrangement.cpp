// Splitting a trimmed face along curves: planar arrangement in the parameter plane.
#include <algorithm>
#include <cmath>
#include <numbers>
#include <optional>

#include "nurbforge/analysis/closest_point.hpp"
#include "nurbforge/analysis/intersect.hpp"
#include "nurbforge/core/error.hpp"
#include "nurbforge/topology/trimmed_surface.hpp"

namespace nf {

namespace {

using P2 = std::array<double, 2>;

P2 at2(const NurbsCurve& c, double t) {
  const Vec3 p = c.evaluate(t);
  return {p.x, p.y};
}

struct Graph {
  std::vector<P2> vertices;
  struct Edge {
    NurbsCurve curve;  // runs from vertex a to vertex b
    std::size_t a, b;
    bool alive{true};
  };
  std::vector<Edge> edges;
  double eps;

  std::size_t vertex(const P2& p) {
    for (std::size_t i = 0; i < vertices.size(); ++i)
      if (std::hypot(vertices[i][0] - p[0], vertices[i][1] - p[1]) <= eps) return i;
    vertices.push_back(p);
    return vertices.size() - 1;
  }
};

// Direction in which the curve leaves its start (forward) or its end (backward), measured by a
// short chord (robust also for tangential departures).
double departure_angle(const NurbsCurve& c, bool from_start) {
  const double a = c.domain_start(), b = c.domain_end(), d = 1e-3 * (b - a);
  const P2 p = at2(c, from_start ? a : b), q = at2(c, from_start ? a + d : b - d);
  return std::atan2(q[1] - p[1], q[0] - p[0]);
}

double polyline_area(const std::vector<P2>& poly) {
  double s = 0.0;
  for (std::size_t i = 0, j = poly.size() - 1; i < poly.size(); j = i++)
    s += (poly[j][0] - poly[i][0]) * (poly[j][1] + poly[i][1]);
  return 0.5 * s;
}

bool inside_polygon(const std::vector<P2>& poly, double u, double v) {
  bool in = false;
  for (std::size_t i = 0, j = poly.size() - 1; i < poly.size(); j = i++) {
    const auto& p = poly[i];
    const auto& q = poly[j];
    if ((p[1] > v) != (q[1] > v) && u < (q[0] - p[0]) * (v - p[1]) / (q[1] - p[1]) + p[0]) in = !in;
  }
  return in;
}

}  // namespace

std::vector<TrimmedSurface> split(const TrimmedSurface& face, const std::vector<NurbsCurve>& cuts,
                                  const Tolerance& tol) {
  const NurbsSurface& surf = face.surface();
  const auto du = surf.domain_u(), dv = surf.domain_v();
  const double size = std::max(du[1] - du[0], dv[1] - dv[0]);
  Tolerance ptol = tol;
  ptol.length = std::max(1e-9 * size, tol.parametric);  // parameter-plane intersections

  // all curves: loop curves and cuts, each with its split parameters
  std::vector<NurbsCurve> curves;
  for (const auto& c : face.outer().curves) curves.push_back(c);
  for (const auto& h : face.holes())
    for (const auto& c : h.curves) curves.push_back(c);
  const std::size_t n_loop = curves.size();
  for (const auto& c : cuts) curves.push_back(c);
  if (cuts.empty()) return {face};

  std::vector<std::vector<double>> params(curves.size());
  for (std::size_t i = 0; i < curves.size(); ++i)
    params[i] = {curves[i].domain_start(), curves[i].domain_end()};
  for (std::size_t i = 0; i < curves.size(); ++i)
    for (std::size_t j = std::max(i + 1, n_loop); j < curves.size(); ++j)  // loop x cut, cut x cut
      for (const auto& x : intersect(curves[i], curves[j], ptol)) {
        params[i].push_back(x.s);
        params[j].push_back(x.t);
      }

  Graph g;
  g.eps = 1e-7 * size;
  for (std::size_t i = 0; i < curves.size(); ++i) {
    auto& ps = params[i];
    std::sort(ps.begin(), ps.end());
    const double span = curves[i].domain_end() - curves[i].domain_start();
    std::vector<double> uniq;
    for (double t : ps)
      if (uniq.empty() || t - uniq.back() > 1e-9 * span) uniq.push_back(t);
    for (std::size_t k = 0; k + 1 < uniq.size(); ++k) {
      NurbsCurve piece = curves[i].sub_curve(uniq[k], uniq[k + 1]);
      const P2 a = at2(piece, piece.domain_start()), b = at2(piece, piece.domain_end());
      if (std::hypot(a[0] - b[0], a[1] - b[1]) <= g.eps) {
        // closed piece (e.g. a full circle): split it once more so that it has two vertices
        const double m = 0.5 * (uniq[k] + uniq[k + 1]);
        for (const auto& [s0, s1] : {std::pair{uniq[k], m}, std::pair{m, uniq[k + 1]}}) {
          NurbsCurve half = curves[i].sub_curve(s0, s1);
          const std::size_t va = g.vertex(at2(half, half.domain_start())), vb = g.vertex(at2(half, half.domain_end()));
          if (i >= n_loop && !face.contains(at2(half, 0.5 * (s0 + s1))[0], at2(half, 0.5 * (s0 + s1))[1])) continue;
          g.edges.push_back({std::move(half), va, vb});
        }
        continue;
      }
      // cut pieces outside the trimmed region are dropped
      if (i >= n_loop) {
        const P2 mid = at2(piece, 0.5 * (piece.domain_start() + piece.domain_end()));
        if (!face.contains(mid[0], mid[1])) continue;
      }
      const std::size_t va = g.vertex(a), vb = g.vertex(b);
      if (va == vb) continue;
      g.edges.push_back({std::move(piece), va, vb});
    }
  }

  // remove dangling edges (vertices of degree 1) iteratively
  for (bool changed = true; changed;) {
    changed = false;
    std::vector<int> degree(g.vertices.size(), 0);
    for (const auto& e : g.edges)
      if (e.alive) {
        ++degree[e.a];
        ++degree[e.b];
      }
    for (auto& e : g.edges)
      if (e.alive && (degree[e.a] < 2 || degree[e.b] < 2)) {
        e.alive = false;
        changed = true;
      }
  }

  // half-edges: 2k forward (a -> b), 2k+1 backward (b -> a)
  const std::size_t H = 2 * g.edges.size();
  const auto origin = [&](std::size_t h) { return h % 2 == 0 ? g.edges[h / 2].a : g.edges[h / 2].b; };
  const auto target = [&](std::size_t h) { return h % 2 == 0 ? g.edges[h / 2].b : g.edges[h / 2].a; };
  const auto angle = [&](std::size_t h) { return departure_angle(g.edges[h / 2].curve, h % 2 == 0); };
  std::vector<std::vector<std::size_t>> out(g.vertices.size());
  for (std::size_t h = 0; h < H; ++h)
    if (g.edges[h / 2].alive) out[origin(h)].push_back(h);
  for (auto& o : out) std::sort(o.begin(), o.end(), [&](std::size_t x, std::size_t y) { return angle(x) < angle(y); });
  // next(h): at target(h), the outgoing half-edge immediately clockwise from twin(h)
  const auto next = [&](std::size_t h) {
    const std::size_t twin = h ^ 1u;
    const auto& o = out[target(h)];
    const auto it = std::find(o.begin(), o.end(), twin);
    const std::size_t k = static_cast<std::size_t>(it - o.begin());
    return o[(k + o.size() - 1) % o.size()];
  };

  struct Cycle {
    std::vector<std::size_t> half_edges;
    std::vector<P2> poly;
    double area;
  };
  std::vector<Cycle> cycles;
  std::vector<bool> used(H, false);
  for (std::size_t h0 = 0; h0 < H; ++h0) {
    if (used[h0] || !g.edges[h0 / 2].alive) continue;
    Cycle c;
    std::size_t h = h0;
    for (std::size_t guard = 0; guard <= H && !used[h]; ++guard) {
      used[h] = true;
      c.half_edges.push_back(h);
      const NurbsCurve& cu = g.edges[h / 2].curve;
      for (int q = 0; q < 8; ++q) {  // polyline in traversal direction
        const double f = q / 8.0;
        const double t = h % 2 == 0 ? cu.domain_start() + f * (cu.domain_end() - cu.domain_start())
                                    : cu.domain_end() - f * (cu.domain_end() - cu.domain_start());
        c.poly.push_back(at2(cu, t));
      }
      h = next(h);
    }
    c.area = polyline_area(c.poly);
    cycles.push_back(std::move(c));
  }

  const auto loop_of = [&](const Cycle& c) {
    TrimLoop l;
    for (std::size_t h : c.half_edges) {
      const NurbsCurve& cu = g.edges[h / 2].curve;
      l.curves.push_back(h % 2 == 0 ? cu : cu.reversed());
    }
    return l;
  };
  // regions: counter-clockwise cycles inside the original face; holes: clockwise cycles assigned
  // to the smallest region containing them
  std::vector<std::size_t> regions;
  for (std::size_t i = 0; i < cycles.size(); ++i) {
    if (cycles[i].area <= 0.0) continue;
    // a point safely inside the cycle: try polyline midpoints nudged to the left
    std::optional<P2> probe;
    const auto& poly = cycles[i].poly;
    for (std::size_t k = 0; k + 1 < poly.size() && !probe; ++k) {
      const double ex = poly[k + 1][0] - poly[k][0], ey = poly[k + 1][1] - poly[k][1];
      const double len = std::hypot(ex, ey);
      if (len <= 0.0) continue;
      const P2 p{0.5 * (poly[k][0] + poly[k + 1][0]) - 1e-4 * size * ey / len,
                 0.5 * (poly[k][1] + poly[k + 1][1]) + 1e-4 * size * ex / len};
      if (inside_polygon(poly, p[0], p[1])) probe = p;
    }
    if (probe && face.contains((*probe)[0], (*probe)[1])) regions.push_back(i);
  }
  std::vector<std::vector<TrimLoop>> holes(regions.size());
  for (std::size_t i = 0; i < cycles.size(); ++i) {
    if (cycles[i].area >= 0.0) continue;
    const P2 p = cycles[i].poly.front();
    std::optional<std::size_t> best;
    for (std::size_t r = 0; r < regions.size(); ++r) {
      const auto& rc = cycles[regions[r]];
      if (-cycles[i].area < rc.area && inside_polygon(rc.poly, p[0], p[1]) &&
          (!best || rc.area < cycles[regions[*best]].area))
        best = r;
    }
    if (best) holes[*best].push_back(loop_of(cycles[i]));
  }
  std::vector<TrimmedSurface> result;
  for (std::size_t r = 0; r < regions.size(); ++r)
    result.emplace_back(surf, loop_of(cycles[regions[r]]), std::move(holes[r]), tol);
  return result;
}

}  // namespace nf
