#include "nurbforge/topology/solid_ops.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

#include "nurbforge/core/error.hpp"
#include "nurbforge/geometry/primitives.hpp"

namespace nf {

namespace {

using P2 = std::array<double, 2>;

// A closed chain of free edges, as 3D curves running head to tail.
struct Opening {
  std::vector<NurbsCurve> curves;
  std::vector<Vec3> samples;  // points along the chain
  Vec3 centroid, normal;      // best-fit plane (Newell normal, unit)
};

std::vector<Opening> free_loops(const Shell& shell, double tol) {
  const auto& E = shell.edges();
  const auto free = shell.free_edges();
  std::vector<bool> used(E.size(), true);
  for (const auto e : free) used[e] = false;
  std::vector<Opening> out;
  for (const auto first : free) {
    if (used[first]) continue;
    Opening op;
    used[first] = true;
    op.curves.push_back(E[first].curve);
    const std::size_t start = E[first].v0;
    std::size_t at = E[first].v1;
    while (at != start) {
      bool found = false;
      for (const auto e : free) {
        if (used[e] || (E[e].v0 != at && E[e].v1 != at)) continue;
        used[e] = true;
        const bool fwd = E[e].v0 == at;
        op.curves.push_back(fwd ? E[e].curve : E[e].curve.reversed());
        at = fwd ? E[e].v1 : E[e].v0;
        found = true;
        break;
      }
      if (!found) throw GeometryError("cap: the free edges do not form closed loops");
    }
    for (const auto& c : op.curves)
      for (int k = 0; k < 16; ++k)
        op.samples.push_back(
            c.evaluate(c.domain_start() + (c.domain_end() - c.domain_start()) * k / 16.0));
    Vec3 cen{}, nrm{};
    for (std::size_t i = 0; i < op.samples.size(); ++i) {
      const Vec3& a = op.samples[i];
      const Vec3& b = op.samples[(i + 1) % op.samples.size()];
      cen = cen + a;
      nrm = nrm +
            Vec3{(a.y - b.y) * (a.z + b.z), (a.z - b.z) * (a.x + b.x), (a.x - b.x) * (a.y + b.y)};
    }
    op.centroid = cen * (1.0 / static_cast<double>(op.samples.size()));
    if (!(norm(nrm) > 0.0)) throw GeometryError("cap: degenerate opening (zero area)");
    op.normal = normalized(nrm);
    double dev = 0.0;
    for (const auto& p : op.samples) dev = std::max(dev, std::abs(dot(p - op.centroid, op.normal)));
    if (dev > tol)
      throw GeometryError("cap: an opening is not planar (deviation " + std::to_string(dev) +
                          " > tolerance)");
    out.push_back(std::move(op));
  }
  return out;
}

double shoelace(const std::vector<P2>& p) {
  double a = 0.0;
  for (std::size_t i = 0, j = p.size() - 1; i < p.size(); j = i++)
    a += p[j][0] * p[i][1] - p[i][0] * p[j][1];
  return 0.5 * a;
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

}  // namespace

Shell cap(const Shell& shell, double tolerance, const std::string& label) {
  if (shell.is_closed()) return shell;
  auto ops = free_loops(shell, tolerance);

  // group openings lying in the same plane
  std::vector<std::vector<std::size_t>> groups;
  for (std::size_t i = 0; i < ops.size(); ++i) {
    bool placed = false;
    for (auto& g : groups) {
      const Opening& r = ops[g.front()];
      if (norm(cross(r.normal, ops[i].normal)) <= 1e-6 &&
          std::abs(dot(ops[i].centroid - r.centroid, r.normal)) <= tolerance) {
        g.push_back(i);
        placed = true;
        break;
      }
    }
    if (!placed) groups.push_back({i});
  }

  std::vector<TrimmedSurface> faces = shell.faces();
  for (const auto& g : groups) {
    const Vec3 n = ops[g.front()].normal;
    const Vec3 s = normalized(any_perpendicular(n)), t = cross(n, s);
    const Vec3 o0 = ops[g.front()].centroid;
    double s0 = 1e300, s1 = -1e300, t0 = 1e300, t1 = -1e300;
    for (const auto i : g)
      for (const auto& c : ops[i].curves)
        for (const auto& p : c.control_points()) {
          s0 = std::min(s0, dot(p - o0, s));
          s1 = std::max(s1, dot(p - o0, s));
          t0 = std::min(t0, dot(p - o0, t));
          t1 = std::max(t1, dot(p - o0, t));
        }
    const double pad = 0.05 * std::max(s1 - s0, t1 - t0);
    s0 -= pad, s1 += pad, t0 -= pad, t1 += pad;
    const Vec3 origin = o0 + s * s0 + t * t0;
    const double W = s1 - s0, H = t1 - t0;
    const NurbsSurface rect = make_rectangle(origin, s * W, t * H);
    const auto to_uv = [&](const Vec3& p) {
      return Vec3{dot(p - origin, s) / W, dot(p - origin, t) / H, 0.0};
    };
    // trim loops: the edge curves mapped affinely (exact); consecutive ends made identical
    std::vector<TrimLoop> loops;
    std::vector<std::vector<P2>> polys;
    for (const auto i : g) {
      TrimLoop loop;
      for (const auto& c3 : ops[i].curves) {
        const NurbsCurve c = c3.is_clamped() ? c3 : c3.clamped();
        std::vector<Vec3> pts;
        for (const auto& p : c.control_points()) pts.push_back(to_uv(p));
        loop.curves.emplace_back(c.degree(), c.knots(), std::move(pts), c.weights());
      }
      for (std::size_t k = 0; k < loop.curves.size(); ++k) {  // close the gaps exactly
        auto& next = loop.curves[(k + 1) % loop.curves.size()];
        auto pts = next.control_points();
        pts.front() = loop.curves[k].control_points().back();
        next = NurbsCurve(next.degree(), next.knots(), std::move(pts), next.weights());
      }
      std::vector<P2> poly;
      for (const auto& p : loop.polyline(16)) poly.push_back(p);
      if (shoelace(poly) < 0.0) {
        loop = loop.reversed();
        std::reverse(poly.begin(), poly.end());
      }
      loops.push_back(std::move(loop));
      polys.push_back(std::move(poly));
    }
    // nesting: depth = number of loops containing a loop; even depth -> outer, odd -> hole
    const std::size_t m = loops.size();
    std::vector<std::size_t> depth(m, 0);
    std::vector<std::ptrdiff_t> parent(m, -1);
    for (std::size_t i = 0; i < m; ++i) {
      double best = 1e300;
      for (std::size_t j = 0; j < m; ++j) {
        if (i == j || !inside(polys[j], polys[i].front())) continue;
        ++depth[i];
        const double a = std::abs(shoelace(polys[j]));
        if (a < best) {
          best = a;
          parent[i] = static_cast<std::ptrdiff_t>(j);
        }
      }
    }
    for (std::size_t i = 0; i < m; ++i) {
      if (depth[i] % 2 != 0) continue;
      std::vector<TrimLoop> holes;
      for (std::size_t j = 0; j < m; ++j)
        if (depth[j] % 2 != 0 && parent[j] == static_cast<std::ptrdiff_t>(i))
          holes.push_back(loops[j].reversed());
      faces.push_back(TrimmedSurface(rect, loops[i], std::move(holes)).with_label(label));
    }
  }
  return Shell::sew(std::move(faces), tolerance).oriented();
}

}  // namespace nf
