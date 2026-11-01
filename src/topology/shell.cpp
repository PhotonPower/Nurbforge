#include "nurbforge/topology/shell.hpp"

#include <algorithm>
#include <cmath>
#include <queue>
#include <string>

#include "nurbforge/analysis/closest_point.hpp"
#include "nurbforge/analysis/properties.hpp"
#include "nurbforge/core/error.hpp"

namespace nf {

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
  for (std::size_t f = 0; f < sh.faces_.size(); ++f) {
    const auto loops = sh.faces_[f].boundary_curves(0.01 * tolerance);
    for (std::size_t l = 0; l < loops.size(); ++l)
      for (std::size_t k = 0; k < loops[l].size(); ++k) {
        const NurbsCurve& c = loops[l][k];
        if (control_diagonal(c) <= tolerance) continue;  // degenerate (pole)
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
  Shell sh = sew(std::move(faces), tol_);
  if (sh.is_closed() && sh.is_consistently_oriented() && sh.volume(1e-6) < 0.0) {
    std::vector<TrimmedSurface> outward;
    for (const auto& t : sh.faces_) outward.push_back(t.flipped());
    sh = sew(std::move(outward), tol_);
  }
  return sh;
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
