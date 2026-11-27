#include "nurbforge/topology/solid_ops.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <numeric>
#include <optional>
#include <string>
#include <vector>

#include "loops.hpp"
#include "nurbforge/analysis/closest_point.hpp"
#include "nurbforge/core/error.hpp"
#include "nurbforge/geometry/primitives.hpp"
#include "nurbforge/modeling/construct.hpp"
#include "nurbforge/modeling/offset.hpp"
#include "nurbforge/modeling/project.hpp"
#include "nurbforge/modeling/surface_tools.hpp"
#include "nurbforge/topology/boolean.hpp"

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

// A planar parallelogram patch S(u, v) = o + (u - u0) / du * e1 + (v - v0) / dv * e2.
struct Parallelogram {
  Vec3 o, e1, e2;
  double u0{0.0}, du{1.0}, v0{0.0}, dv{1.0};
};

std::optional<Parallelogram> parallelogram(const NurbsSurface& s, double tol) {
  if (s.degree_u() != 1 || s.degree_v() != 1 || s.num_u() != 2 || s.num_v() != 2) return {};
  const auto& w = s.weights();
  for (const double x : w)
    if (std::abs(x - w.front()) > 1e-12 * std::abs(w.front())) return {};
  const auto ku = s.knots_u(), kv = s.knots_v();
  Parallelogram p;
  p.o = s.control_point(0, 0);
  p.e1 = s.control_point(1, 0) - p.o;
  p.e2 = s.control_point(0, 1) - p.o;
  if (distance(s.control_point(1, 1), p.o + p.e1 + p.e2) > tol) return {};
  p.u0 = ku.front();
  p.du = ku.back() - ku.front();
  p.v0 = kv.front();
  p.dv = kv.back() - kv.front();
  return p;
}

bool same_surface(const NurbsSurface& a, const NurbsSurface& b) {
  if (a.degree_u() != b.degree_u() || a.degree_v() != b.degree_v() || a.num_u() != b.num_u() ||
      a.num_v() != b.num_v() || a.knots_u() != b.knots_u() || a.knots_v() != b.knots_v() ||
      a.weights() != b.weights())
    return false;
  const auto &P = a.control_points(), &Q = b.control_points();
  for (std::size_t k = 0; k < P.size(); ++k)
    if (P[k].x != Q[k].x || P[k].y != Q[k].y || P[k].z != Q[k].z) return false;
  return true;
}

NurbsCurve map_points(const NurbsCurve& c, const std::function<Vec3(const Vec3&)>& f) {
  std::vector<Vec3> pts;
  for (const auto& p : c.control_points()) pts.push_back(f(p));
  return {c.degree(), c.knots(), std::move(pts), c.weights()};
}

Vec3 first(const NurbsCurve& c) { return c.evaluate(c.domain_start()); }
Vec3 last(const NurbsCurve& c) { return c.evaluate(c.domain_end()); }

}  // namespace

namespace detail {

// The region bounded by curves with the region on their left (outer loops counter-clockwise,
// holes clockwise), all in one parameter plane: shared pieces (opposite directions) cancel, the
// rest is chained into loops. Empty if the curves do not form closed loops.
std::vector<TrimLoop> union_loops(std::vector<NurbsCurve> cs, double eps) {
  // T-junctions: split every curve at the end points of other curves lying inside it
  std::vector<NurbsCurve> pieces;
  for (std::size_t i = 0; i < cs.size(); ++i) {
    const NurbsCurve& c = cs[i];
    const double a = c.domain_start(), b = c.domain_end();
    std::vector<double> cut{a, b};
    for (std::size_t j = 0; j < cs.size(); ++j) {
      if (j == i) continue;
      for (const Vec3& p : {first(cs[j]), last(cs[j])}) {
        const CurvePoint cp = closest_point(c, p);
        if (cp.distance <= eps && cp.u > a + 1e-9 * (b - a) && cp.u < b - 1e-9 * (b - a))
          cut.push_back(cp.u);
      }
    }
    std::sort(cut.begin(), cut.end());
    for (std::size_t k = 0; k + 1 < cut.size(); ++k)
      if (cut[k + 1] - cut[k] > 1e-9 * (b - a))
        pieces.push_back(cut.size() == 2 ? c : c.sub_curve(cut[k], cut[k + 1]));
  }
  // cancel pairs running along each other in opposite directions
  std::vector<bool> gone(pieces.size(), false);
  for (std::size_t i = 0; i < pieces.size(); ++i)
    for (std::size_t j = i + 1; j < pieces.size() && !gone[i]; ++j) {
      if (gone[j] || distance(first(pieces[i]), last(pieces[j])) > eps ||
          distance(last(pieces[i]), first(pieces[j])) > eps)
        continue;
      const NurbsCurve& c = pieces[i];
      const Vec3 mid = c.evaluate(0.5 * (c.domain_start() + c.domain_end()));
      if (closest_point(pieces[j], mid).distance <= eps) gone[i] = gone[j] = true;
    }
  std::vector<NurbsCurve> rest;
  for (std::size_t i = 0; i < pieces.size(); ++i)
    if (!gone[i]) rest.push_back(std::move(pieces[i]));
  // chain into loops
  std::vector<TrimLoop> loops;
  std::vector<bool> used(rest.size(), false);
  for (std::size_t s0 = 0; s0 < rest.size(); ++s0) {
    if (used[s0]) continue;
    TrimLoop loop;
    used[s0] = true;
    loop.curves.push_back(rest[s0]);
    while (distance(last(loop.curves.back()), first(loop.curves.front())) > eps) {
      std::size_t next = rest.size();
      for (std::size_t j = 0; j < rest.size(); ++j)
        if (!used[j] && distance(first(rest[j]), last(loop.curves.back())) <= eps) {
          next = j;
          break;
        }
      if (next == rest.size()) return {};  // open chain
      used[next] = true;
      loop.curves.push_back(rest[next]);
    }
    // make consecutive ends identical
    for (std::size_t k = 0; k < loop.curves.size(); ++k) {
      auto& nx = loop.curves[(k + 1) % loop.curves.size()];
      auto pts = nx.control_points();
      pts.front() = loop.curves[k].control_points().back();
      nx = NurbsCurve(nx.degree(), nx.knots(), std::move(pts), nx.weights());
    }
    loops.push_back(std::move(loop));
  }
  return loops;
}

}  // namespace detail

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

Shell thicken(const Shell& shell, double thickness, bool symmetric, double tolerance,
              const std::string& wall_label) {
  if (!(std::abs(thickness) > 0.0)) throw GeometryError("thicken: thickness must not be zero");
  if (shell.is_closed()) throw GeometryError("thicken: the shell is closed (nothing to thicken)");
  const double d0 = symmetric ? -0.5 * thickness : 0.0,
               d1 = symmetric ? 0.5 * thickness : thickness;
  const double fit = 0.01 * tolerance;
  const auto offset_face = [&](const TrimmedSurface& f, double d, const std::string& suffix) {
    if (d == 0.0) return f;
    const NurbsSurface s = offset_surface(f.surface(), d, fit);
    return TrimmedSurface(s, f.outer(), f.holes())
        .with_label(f.label().empty() ? suffix.substr(1) : f.label() + suffix);
  };
  // 3D pieces of the free edges, to recognise the loop curves lying on them
  std::vector<const NurbsCurve*> free;
  for (const auto e : shell.free_edges()) free.push_back(&shell.edges()[e].curve);
  const auto on_free_edge = [&](const NurbsCurve& c) {
    for (const double s : {0.25, 0.5, 0.75}) {
      const Vec3 p = c.evaluate(c.domain_start() + s * (c.domain_end() - c.domain_start()));
      bool near = false;
      for (const NurbsCurve* e : free)
        if (closest_point(*e, p).distance <= tolerance) {
          near = true;
          break;
        }
      if (!near) return false;
    }
    return true;
  };
  std::vector<TrimmedSurface> faces;
  for (const auto& f : shell.faces()) {
    const TrimmedSurface a = offset_face(f, d0, symmetric ? "_offset" : "");
    const TrimmedSurface b = offset_face(f, d1, "_offset");
    const auto ca = a.boundary_curves(fit), cb = b.boundary_curves(fit),
               c0 = f.boundary_curves(fit);
    for (std::size_t l = 0; l < c0.size(); ++l)
      for (std::size_t k = 0; k < c0[l].size(); ++k)
        if (on_free_edge(c0[l][k]))
          faces.push_back(TrimmedSurface(ruled(ca[l][k], cb[l][k])).with_label(wall_label));
    faces.push_back(a);
    faces.push_back(b);
  }
  Shell out = Shell::sew(std::move(faces), tolerance);
  if (!out.is_closed())
    throw GeometryError("thicken: the result is not closed (" +
                        std::to_string(out.free_edges().size()) +
                        " free edges) - the faces must join smoothly so that their offsets meet");
  return out.oriented();
}

Shell make_pipe(const NurbsCurve& path, const std::function<double(double)>& radius,
                double tolerance) {
  SweepOptions opt;
  opt.tolerance = tolerance;
  const NurbsSurface wall = pipe_surface(path, radius, opt);
  const double sew = 100.0 * tolerance;
  const Vec3 p0 = path.evaluate(path.domain_start()), p1 = path.evaluate(path.domain_end());
  Shell s = Shell::sew({TrimmedSurface(wall).with_label("wall")}, sew);
  if (distance(p0, p1) > sew) {
    s = cap(s, sew, "end");
    // name the end faces by the path end they close
    std::vector<TrimmedSurface> faces;
    for (const auto& f : s.faces()) {
      if (f.label() != "end") {
        faces.push_back(f);
        continue;
      }
      const auto c = f.boundary_curves(0.01 * tolerance).front().front();
      const Vec3 q = c.evaluate(c.domain_start());
      faces.push_back(f.with_label(distance(q, p0) <= distance(q, p1) ? "inlet" : "outlet"));
    }
    s = Shell::sew(std::move(faces), sew);
  }
  return s.oriented();
}

Shell make_pipe(const NurbsCurve& path, double start_radius, double end_radius, double tolerance) {
  const double r1 = end_radius < 0.0 ? start_radius : end_radius;
  return make_pipe(
      path, [=](double t) { return start_radius + (r1 - start_radius) * t; }, tolerance);
}

Shell remove_holes(const Shell& shell, double max_size, double tolerance) {
  if (!shell.is_closed()) throw GeometryError("remove_holes: the shell must be closed");
  const auto extent = [&](const std::vector<NurbsCurve>& loop) {
    Vec3 lo{1e300, 1e300, 1e300}, hi{-1e300, -1e300, -1e300};
    for (const auto& c : loop)
      for (const auto& p : c.control_points()) {
        lo = {std::min(lo.x, p.x), std::min(lo.y, p.y), std::min(lo.z, p.z)};
        hi = {std::max(hi.x, p.x), std::max(hi.y, p.y), std::max(hi.z, p.z)};
      }
    return distance(lo, hi);
  };
  std::vector<TrimmedSurface> faces;
  std::size_t dropped = 0;
  for (const auto& f : shell.faces()) {
    const auto loops = f.boundary_curves(0.01 * tolerance);
    std::vector<TrimLoop> keep;
    for (std::size_t h = 0; h < f.holes().size(); ++h) {
      if (extent(loops[h + 1]) < max_size)
        ++dropped;
      else
        keep.push_back(f.holes()[h]);
    }
    faces.push_back(
        keep.size() == f.holes().size()
            ? f
            : TrimmedSurface(f.surface(), f.outer(), std::move(keep)).with_label(f.label()));
  }
  if (dropped == 0) return shell;
  // remove the faces of the features: faces with free edges, repeatedly
  for (std::size_t round = 0; round < shell.faces().size(); ++round) {
    const Shell s = Shell::sew(faces, tolerance);
    const auto free = s.free_edges();
    if (free.empty()) return s.oriented();
    std::vector<bool> drop(faces.size(), false);
    for (const auto e : free)
      for (const auto& use : s.edges()[e].uses) drop[use.face] = true;
    std::vector<TrimmedSurface> rest;
    for (std::size_t k = 0; k < faces.size(); ++k)
      if (!drop[k]) rest.push_back(s.faces()[k]);
    if (rest.size() == faces.size() || rest.empty()) break;
    faces = std::move(rest);
  }
  throw GeometryError(
      "remove_holes: the shell does not close after removing the features "
      "(only holes - through or blind - are supported)");
}

Shell hollow(const Shell& shell, double thickness, const std::vector<std::string>& open_labels,
             double tolerance) {
  if (!(thickness > 0.0)) throw GeometryError("hollow: thickness must be positive");
  if (!shell.is_closed()) throw GeometryError("hollow: the shell must be closed");
  const Shell outer = shell.oriented();
  if (!open_labels.empty()) {
    std::vector<TrimmedSurface> rest;
    for (const auto& f : outer.faces())
      if (std::find(open_labels.begin(), open_labels.end(), f.label()) == open_labels.end())
        rest.push_back(f);
    if (rest.size() == outer.faces().size())
      throw GeometryError("hollow: no face carries one of the open labels");
    // outward normals: thicken against them (inwards)
    return thicken(Shell::sew(std::move(rest), tolerance), -thickness, false, tolerance);
  }
  std::vector<TrimmedSurface> faces = outer.faces();
  for (const auto& f : outer.faces()) {
    const NurbsSurface s = offset_surface(f.surface(), -thickness, 0.01 * tolerance);
    faces.push_back(TrimmedSurface(s, f.outer(), f.holes())
                        .with_label(f.label().empty() ? "inner" : f.label() + "_inner"));
  }
  const Shell out = Shell::sew(std::move(faces), tolerance);
  if (!out.is_closed() || out.components().size() != 2)
    throw GeometryError(
        "hollow: the offset faces do not form a closed inner shell - the faces "
        "must join smoothly so that their offsets meet");
  return out.oriented();
}

Shell merge_faces(const Shell& shell, double tolerance) {
  const auto& F = shell.faces();
  const std::size_t n = F.size();
  std::vector<std::optional<Parallelogram>> plane(n);
  for (std::size_t i = 0; i < n; ++i) plane[i] = parallelogram(F[i].surface(), tolerance);
  const auto compatible = [&](std::size_t i, std::size_t j) {
    if (F[i].label() != F[j].label()) return false;
    if (same_surface(F[i].surface(), F[j].surface())) return true;
    if (!plane[i] || !plane[j]) return false;
    const Vec3 ni = normalized(cross(plane[i]->e1, plane[i]->e2));
    const Vec3 nj = normalized(cross(plane[j]->e1, plane[j]->e2));
    return dot(ni, nj) > 1.0 - 1e-9 && std::abs(dot(plane[j]->o - plane[i]->o, ni)) <= tolerance;
  };
  // groups of compatible faces connected through shared edges (union-find)
  std::vector<std::size_t> root(n);
  std::iota(root.begin(), root.end(), 0);
  const auto find = [&](std::size_t x) {
    while (root[x] != x) x = root[x] = root[root[x]];
    return x;
  };
  for (const auto& e : shell.edges())
    if (e.uses.size() == 2 && e.uses[0].face != e.uses[1].face &&
        compatible(e.uses[0].face, e.uses[1].face))
      root[find(e.uses[0].face)] = find(e.uses[1].face);
  std::vector<std::vector<std::size_t>> groups(n);
  for (std::size_t i = 0; i < n; ++i) groups[find(i)].push_back(i);

  std::vector<TrimmedSurface> out;
  for (const auto& g : groups) {
    if (g.empty()) continue;
    if (g.size() == 1) {
      out.push_back(F[g.front()]);
      continue;
    }
    // carrier surface and the map of every member's trim curves into its parameter plane
    NurbsSurface carrier = F[g.front()].surface();
    std::function<Vec3(std::size_t, const Vec3&)> to_carrier = [](std::size_t, const Vec3& p) {
      return p;
    };
    double eps = 1e-9;
    bool all_same = true;
    for (const auto i : g) all_same = all_same && same_surface(F[i].surface(), carrier);
    if (all_same) {
      const auto du = carrier.domain_u(), dv = carrier.domain_v();
      eps = 1e-8 * std::max(du[1] - du[0], dv[1] - dv[0]);
    } else {
      const Parallelogram& p0 = *plane[g.front()];
      const Vec3 nrm = normalized(cross(p0.e1, p0.e2));
      const Vec3 sx = normalized(p0.e1), ty = cross(nrm, sx);
      const auto world = [&](std::size_t i, const Vec3& q) {
        const Parallelogram& p = *plane[i];
        return p.o + ((q.x - p.u0) / p.du) * p.e1 + ((q.y - p.v0) / p.dv) * p.e2;
      };
      double s0 = 1e300, s1 = -1e300, t0 = 1e300, t1 = -1e300;
      for (const auto i : g)
        for (const auto& c : F[i].outer().curves)
          for (const auto& q : c.control_points()) {
            const Vec3 x = world(i, q) - p0.o;
            s0 = std::min(s0, dot(x, sx)), s1 = std::max(s1, dot(x, sx));
            t0 = std::min(t0, dot(x, ty)), t1 = std::max(t1, dot(x, ty));
          }
      const double W = s1 - s0, H = t1 - t0;
      const Vec3 origin = p0.o + sx * s0 + ty * t0;
      carrier = make_rectangle(origin, sx * W, ty * H);
      to_carrier = [=](std::size_t i, const Vec3& q) {
        const Vec3 x = world(i, q) - origin;
        return Vec3{dot(x, sx) / W, dot(x, ty) / H, 0.0};
      };
      eps = tolerance / std::max(W, H);
    }
    std::vector<NurbsCurve> curves;
    for (const auto i : g) {
      const auto map = [&](const NurbsCurve& c) {
        return map_points(c, [&](const Vec3& q) { return to_carrier(i, q); });
      };
      for (const auto& c : F[i].outer().curves) curves.push_back(map(c));
      for (const auto& h : F[i].holes())
        for (const auto& c : h.curves) curves.push_back(map(c));
    }
    std::vector<TrimmedSurface> merged;
    try {
      auto loops = detail::union_loops(std::move(curves), eps);
      std::vector<TrimLoop> outers, holes;
      for (auto& l : loops) (l.signed_area() > 0.0 ? outers : holes).push_back(std::move(l));
      if (outers.empty()) throw GeometryError("no outer loop");
      std::vector<std::vector<TrimLoop>> inner(outers.size());
      for (auto& h : holes) {
        const auto q = h.polyline(4).front();
        std::size_t best = outers.size();
        double area = 1e300;
        for (std::size_t k = 0; k < outers.size(); ++k) {
          const auto poly = outers[k].polyline(8);
          bool in = false;
          for (std::size_t a = 0, b = poly.size() - 1; a < poly.size(); b = a++)
            if ((poly[a][1] > q[1]) != (poly[b][1] > q[1]) &&
                q[0] < poly[b][0] + (poly[a][0] - poly[b][0]) * (q[1] - poly[b][1]) /
                                        (poly[a][1] - poly[b][1]))
              in = !in;
          if (in && outers[k].signed_area() < area) {
            area = outers[k].signed_area();
            best = k;
          }
        }
        if (best == outers.size()) throw GeometryError("hole outside every outer loop");
        inner[best].push_back(std::move(h));
      }
      for (std::size_t k = 0; k < outers.size(); ++k)
        merged.push_back(TrimmedSurface(carrier, outers[k], std::move(inner[k]))
                             .with_label(F[g.front()].label()));
    } catch (const GeometryError&) {
      merged.clear();
      for (const auto i : g) merged.push_back(F[i]);  // leave the group as it is
    }
    for (auto& m : merged) out.push_back(std::move(m));
  }
  return Shell::sew(std::move(out), tolerance);
}

namespace {

Vec3 unit_direction(const Vec3& d, const char* who) {
  const double n = norm(d);
  if (!(n > 0.0) || !std::isfinite(n))
    throw GeometryError(std::string(who) + ": the direction must be a non-zero vector");
  return d / n;
}

Shell labelled(const Shell& sh, const std::string& label) {
  return sh.with_face_labels(std::vector<std::string>(sh.faces().size(), label));
}

}  // namespace

Shell make_slab(const NurbsCurve& curve, double width, const Vec3& direction, double tolerance) {
  const Vec3 n = unit_direction(direction, "make_slab");
  if (width == 0.0) throw GeometryError("make_slab: width must be non-zero");
  // strips between the pieces of the curve (v = 0) and of its offset (v = 1)
  const std::vector<NurbsSurface> strips = ribbon(curve, width, n, 0.1 * tolerance);
  const bool closed = distance(curve.evaluate(curve.domain_start()),
                               curve.evaluate(curve.domain_end())) <= tolerance;
  std::vector<TrimmedSurface> walls;
  for (std::size_t i = 0; i < strips.size(); ++i) {
    const NurbsSurface& s = strips[i];
    const auto du = s.domain_u(), dv = s.domain_v();
    walls.emplace_back(extrude(s.iso_curve_v(dv[0]), direction));
    walls.emplace_back(extrude(s.iso_curve_v(dv[1]), direction));
    if (!closed && i == 0) walls.emplace_back(extrude(s.iso_curve_u(du[0]), direction));
    if (!closed && i + 1 == strips.size())
      walls.emplace_back(extrude(s.iso_curve_u(du[1]), direction));
  }
  return cap(labelled(Shell::sew(std::move(walls), tolerance), "side"), tolerance);
}

std::vector<Shell> wire_cut(const Shell& shell, const NurbsCurve& curve, const Vec3& direction,
                            double tolerance, const std::string& cut_label) {
  const Vec3 n = unit_direction(direction, "wire_cut");
  if (!shell.is_closed()) throw GeometryError("wire_cut: the shell must be closed");
  // extent of the solid along the direction and overall
  Vec3 lo{1e300, 1e300, 1e300}, hi{-1e300, -1e300, -1e300};
  double h0 = 1e300, h1 = -1e300;
  for (const auto& f : shell.faces())
    for (const auto& p : f.surface().control_points()) {
      lo = {std::min(lo.x, p.x), std::min(lo.y, p.y), std::min(lo.z, p.z)};
      hi = {std::max(hi.x, p.x), std::max(hi.y, p.y), std::max(hi.z, p.z)};
      h0 = std::min(h0, dot(p, n));
      h1 = std::max(h1, dot(p, n));
    }
  const double size = distance(lo, hi);
  const double margin = 0.25 * size + 10.0 * tolerance;
  // the curve in the plane below the solid
  const NurbsCurve base = project_to_plane(curve, Plane{(h0 - margin) * n, n}, n);
  const Vec3 A = base.evaluate(base.domain_start()), B = base.evaluate(base.domain_end());
  std::vector<NurbsCurve> loop{base};
  if (distance(A, B) > tolerance) {
    // closed on the left, beyond the solid and the curve
    const Vec3 left = normalized(cross(n, B - A));
    double reach = 0.0;
    for (const auto& p : base.control_points()) reach = std::max(reach, std::abs(dot(p - A, left)));
    const double w = 2.0 * (size + reach + distance(A, 0.5 * (lo + hi))) + margin;
    loop.push_back(make_line(B, B + w * left));
    loop.push_back(make_line(B + w * left, A + w * left));
    loop.push_back(make_line(A + w * left, A));
  }
  const Vec3 height = (h1 - h0 + 2.0 * margin) * n;
  std::vector<TrimmedSurface> walls;
  for (const auto& c : loop) walls.emplace_back(extrude(c, height));
  const Shell prism =
      cap(labelled(Shell::sew(std::move(walls), tolerance), cut_label), tolerance, cut_label);
  std::vector<Shell> out;
  for (const BooleanOp op : {BooleanOp::Intersection, BooleanOp::Difference}) {
    Shell part = boolean(shell, prism, op, tolerance);
    if (!part.faces().empty()) out.push_back(std::move(part));
  }
  return out;
}

}  // namespace nf
