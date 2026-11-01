#include "nurbforge/topology/shell.hpp"

#include <algorithm>
#include <cmath>
#include <queue>
#include <string>

#include "nurbforge/analysis/closest_point.hpp"
#include "nurbforge/analysis/properties.hpp"
#include "nurbforge/core/error.hpp"
#include "nurbforge/topology/boolean.hpp"

namespace nf {

namespace {

// A point on the shell: midpoint of its first edge, else the centre of the first face's domain.
Vec3 probe_point(const Shell& s) {
  if (!s.edges().empty()) {
    const NurbsCurve& c = s.edges().front().curve;
    return c.evaluate(0.5 * (c.domain_start() + c.domain_end()));
  }
  const NurbsSurface& f = s.faces().front().surface();
  const auto du = f.domain_u(), dv = f.domain_v();
  return f.evaluate(0.5 * (du[0] + du[1]), 0.5 * (dv[0] + dv[1]));
}

}  // namespace

namespace {

double control_diagonal(const NurbsCurve& c) {
  const auto& P = c.control_points();
  Vec3 lo = P.front(), hi = P.front();
  for (const auto& p : P) {
    lo = {std::min(lo.x, p.x), std::min(lo.y, p.y), std::min(lo.z, p.z)};
    hi = {std::max(hi.x, p.x), std::max(hi.y, p.y), std::max(hi.z, p.z)};
  }
  return distance(lo, hi);
}

// Every sample of a lies on b within tol.
bool lies_on(const NurbsCurve& a, const NurbsCurve& b, double tol) {
  for (int i = 0; i <= 6; ++i) {
    const double t = a.domain_start() + (a.domain_end() - a.domain_start()) * i / 6.0;
    if (closest_point(b, a.evaluate(t)).distance > tol) return false;
  }
  return true;
}

// Direction of a relative to b: true if opposite (compared at a's midpoint).
bool opposite(const NurbsCurve& a, const NurbsCurve& b) {
  const double tm = 0.5 * (a.domain_start() + a.domain_end());
  const Vec3 ta = a.derivatives(tm, 1)[1];
  const CurvePoint cp = closest_point(b, a.evaluate(tm));
  const Vec3 tb = b.derivatives(cp.u, 1)[1];
  return dot(ta, tb) < 0.0;
}

}  // namespace

Shell Shell::sew(std::vector<TrimmedSurface> faces, double tolerance) {
  Shell sh;
  sh.faces_ = std::move(faces);
  sh.tol_ = tolerance;
  const auto vertex = [&](const Vec3& p) {
    for (std::size_t i = 0; i < sh.vertices_.size(); ++i)
      if (distance(sh.vertices_[i], p) <= tolerance) return i;
    sh.vertices_.push_back(p);
    return sh.vertices_.size() - 1;
  };
  // all non-degenerate boundary curves in 3D
  struct Piece {
    NurbsCurve c;
    std::size_t f, l, k;
  };
  std::vector<Piece> raw;
  for (std::size_t f = 0; f < sh.faces_.size(); ++f) {
    const auto loops = sh.faces_[f].boundary_curves(0.01 * tolerance);
    for (std::size_t l = 0; l < loops.size(); ++l)
      for (std::size_t k = 0; k < loops[l].size(); ++k)
        if (control_diagonal(loops[l][k]) > tolerance) raw.push_back({loops[l][k], f, l, k});
  }
  // resolve T-junctions: split every curve at the end points of other curves lying on it (e.g.
  // a circle split at different points on the two faces sharing it)
  std::vector<Piece> pieces;
  for (std::size_t i = 0; i < raw.size(); ++i) {
    const NurbsCurve& c = raw[i].c;
    const double a = c.domain_start(), b = c.domain_end();
    std::vector<double> cuts;
    for (std::size_t j = 0; j < raw.size(); ++j) {
      if (j == i) continue;
      for (const Vec3& q : {raw[j].c.evaluate(raw[j].c.domain_start()), raw[j].c.evaluate(raw[j].c.domain_end())}) {
        const CurvePoint cp = closest_point(c, q);
        const double margin = 1e-9 * (b - a);
        if (cp.distance <= tolerance && cp.u > a + margin && cp.u < b - margin &&
            distance(cp.point, c.evaluate(a)) > tolerance && distance(cp.point, c.evaluate(b)) > tolerance)
          cuts.push_back(cp.u);
      }
    }
    std::sort(cuts.begin(), cuts.end());
    double from = a;
    for (double t : cuts) {
      if (t - from <= 1e-9 * (b - a)) continue;
      pieces.push_back({c.sub_curve(from, t), raw[i].f, raw[i].l, raw[i].k});
      from = t;
    }
    pieces.push_back({from == a ? c : c.sub_curve(from, b), raw[i].f, raw[i].l, raw[i].k});
  }
  for (const auto& pc : pieces) {
    {
      const NurbsCurve& c = pc.c;
      const std::size_t f = pc.f, l = pc.l, k = pc.k;
      {
        bool found = false;
        for (auto& e : sh.edges_)
          if (lies_on(c, e.curve, tolerance) && lies_on(e.curve, c, tolerance)) {
            e.uses.push_back({f, l, k, opposite(c, e.curve)});
            found = true;
            break;
          }
        if (found) continue;
        ShellEdge e{c, vertex(c.evaluate(c.domain_start())), vertex(c.evaluate(c.domain_end())), {}};
        e.uses.push_back({f, l, k, false});
        sh.edges_.push_back(std::move(e));
      }
    }
  }
  return sh;
}

std::vector<std::size_t> Shell::free_edges() const {
  std::vector<std::size_t> out;
  for (std::size_t i = 0; i < edges_.size(); ++i)
    if (edges_[i].uses.size() == 1) out.push_back(i);
  return out;
}

std::vector<std::size_t> Shell::non_manifold_edges() const {
  std::vector<std::size_t> out;
  for (std::size_t i = 0; i < edges_.size(); ++i)
    if (edges_[i].uses.size() > 2) out.push_back(i);
  return out;
}

bool Shell::is_closed() const {
  return std::all_of(edges_.begin(), edges_.end(), [](const ShellEdge& e) { return e.uses.size() == 2; });
}

bool Shell::is_consistently_oriented() const {
  return std::all_of(edges_.begin(), edges_.end(), [](const ShellEdge& e) {
    return e.uses.size() != 2 || e.uses[0].reversed != e.uses[1].reversed;
  });
}

std::vector<Shell> Shell::components() const {
  std::vector<Shell> out;
  for (const auto& comp : face_components()) {
    std::vector<TrimmedSurface> faces;
    for (std::size_t f : comp) faces.push_back(faces_[f]);
    out.push_back(sew(std::move(faces), tol_));
  }
  return out;
}

std::vector<std::vector<std::size_t>> Shell::face_components() const {
  const std::size_t F = faces_.size();
  std::vector<std::size_t> parent(F);
  for (std::size_t f = 0; f < F; ++f) parent[f] = f;
  const auto find = [&](std::size_t f) {
    while (parent[f] != f) f = parent[f] = parent[parent[f]];
    return f;
  };
  for (const auto& e : edges_)
    for (std::size_t k = 1; k < e.uses.size(); ++k) parent[find(e.uses[k].face)] = find(e.uses[0].face);
  std::vector<std::vector<std::size_t>> comps;
  std::vector<std::size_t> index(F, F);
  for (std::size_t f = 0; f < F; ++f) {
    const std::size_t r = find(f);
    if (index[r] == F) {
      index[r] = comps.size();
      comps.emplace_back();
    }
    comps[index[r]].push_back(f);
  }
  return comps;
}

Shell Shell::oriented() const {
  const std::size_t F = faces_.size();
  std::vector<int> flip(F, -1);  // -1: unvisited, 0/1: keep/flip
  for (std::size_t start = 0; start < F; ++start) {
    if (flip[start] >= 0) continue;
    flip[start] = 0;
    std::queue<std::size_t> todo;
    todo.push(start);
    while (!todo.empty()) {
      const std::size_t f = todo.front();
      todo.pop();
      for (const auto& e : edges_) {
        if (e.uses.size() != 2 || e.uses[0].face == e.uses[1].face) continue;
        for (int side = 0; side < 2; ++side) {
          const Coedge& a = e.uses[static_cast<std::size_t>(side)];
          const Coedge& b = e.uses[static_cast<std::size_t>(1 - side)];
          if (a.face != f) continue;
          // consistent if the two uses run in opposite directions after flipping
          const int need = (a.reversed == b.reversed) ? 1 - flip[f] : flip[f];
          if (flip[b.face] < 0) {
            flip[b.face] = need;
            todo.push(b.face);
          } else if (flip[b.face] != need) {
            throw GeometryError("shell: the faces cannot be oriented consistently (non-orientable)");
          }
        }
      }
    }
  }
  std::vector<TrimmedSurface> faces;
  for (std::size_t f = 0; f < F; ++f) faces.push_back(flip[f] ? faces_[f].flipped() : faces_[f]);
  const auto comps = face_components();
  if (comps.size() == 1) {  // a single component: outward if closed
    Shell sh = sew(std::move(faces), tol_);
    if (sh.is_closed() && sh.is_consistently_oriented() && sh.volume(1e-6) < 0.0) {
      std::vector<TrimmedSurface> outward;
      for (const auto& t : sh.faces_) outward.push_back(t.flipped());
      sh = sew(std::move(outward), tol_);
    }
    return sh;
  }
  // closed components: orientation by nesting depth (even: outer boundary, odd: void)
  std::vector<Shell> parts;
  std::vector<bool> closed(comps.size(), false);
  for (std::size_t c = 0; c < comps.size(); ++c) {
    std::vector<TrimmedSurface> fc;
    for (std::size_t f : comps[c]) fc.push_back(faces[f]);
    parts.push_back(sew(std::move(fc), tol_));
    closed[c] = parts[c].is_closed() && parts[c].is_consistently_oriented();
  }
  for (std::size_t c = 0; c < comps.size(); ++c) {
    if (!closed[c]) continue;
    int depth = 0;
    if (comps.size() > 1) {
      const Vec3 probe = probe_point(parts[c]);
      for (std::size_t d = 0; d < comps.size(); ++d)
        if (d != c && closed[d] && point_in_shell(parts[d], probe, tol_)) ++depth;
    }
    const bool outward = depth % 2 == 0;
    if ((parts[c].volume(1e-6) > 0.0) != outward)
      for (std::size_t f : comps[c]) faces[f] = faces[f].flipped();
  }
  return sew(std::move(faces), tol_);
}

double Shell::area(double rel_tol) const {
  double a = 0.0;
  for (const auto& f : faces_) a += f.area(rel_tol);
  return a;
}

double Shell::volume(double rel_tol) const {
  if (!is_closed()) throw GeometryError("shell: volume requires a closed shell");
  if (!is_consistently_oriented()) throw GeometryError("shell: volume requires a consistently oriented shell");
  // scale for the absolute error target: area times model size
  double size = 0.0;
  for (const auto& v : vertices_) size = std::max(size, norm(v));
  for (const auto& e : edges_) size = std::max(size, control_diagonal(e.curve));
  const double abs_tol = rel_tol * std::max(area(1e-6), 1e-300) * std::max(size, 1e-300);
  double V = 0.0;
  for (const auto& f : faces_) {
    const NurbsSurface& s = f.surface();
    V += f.integrate(
        [&](double u, double v) {
          const auto d = s.first_derivatives(u, v);
          return dot(d[0], cross(d[1], d[2])) / 3.0;
        },
        rel_tol, abs_tol / static_cast<double>(faces_.size()));
  }
  return V;
}

}  // namespace nf
