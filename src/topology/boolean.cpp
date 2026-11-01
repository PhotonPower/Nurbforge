#include "nurbforge/topology/boolean.hpp"

#include <algorithm>
#include <cmath>
#include <optional>

#include "nurbforge/analysis/intersect.hpp"
#include "nurbforge/analysis/surface_intersect.hpp"
#include "nurbforge/core/error.hpp"
#include "nurbforge/geometry/primitives.hpp"

namespace nf {

namespace {

struct Box3 {
  Vec3 lo{1e300, 1e300, 1e300}, hi{-1e300, -1e300, -1e300};
  void add(const Vec3& p) {
    lo = {std::min(lo.x, p.x), std::min(lo.y, p.y), std::min(lo.z, p.z)};
    hi = {std::max(hi.x, p.x), std::max(hi.y, p.y), std::max(hi.z, p.z)};
  }
  [[nodiscard]] bool overlaps(const Box3& o, double eps) const {
    return lo.x <= o.hi.x + eps && o.lo.x <= hi.x + eps && lo.y <= o.hi.y + eps && o.lo.y <= hi.y + eps &&
           lo.z <= o.hi.z + eps && o.lo.z <= hi.z + eps;
  }
};

Box3 face_box(const TrimmedSurface& f) {
  Box3 b;
  for (const auto& p : f.surface().control_points()) b.add(p);
  return b;
}

Box3 shell_box(const Shell& s) {
  Box3 b;
  for (const auto& f : s.faces())
    for (const auto& p : f.surface().control_points()) b.add(p);
  return b;
}

// A parameter point inside the trimmed face (nudged inwards from the outer loop).
std::optional<std::array<double, 2>> interior_point(const TrimmedSurface& f) {
  const auto poly = f.outer().polyline(8);
  const auto du = f.surface().domain_u(), dv = f.surface().domain_v();
  const double size = std::max(du[1] - du[0], dv[1] - dv[0]);
  for (double nudge : {1e-3, 1e-2, 1e-4})
    for (std::size_t k = 0; k < poly.size(); ++k) {
      const auto& p = poly[k];
      const auto& q = poly[(k + 1) % poly.size()];
      const double ex = q[0] - p[0], ey = q[1] - p[1], len = std::hypot(ex, ey);
      if (len <= 0.0) continue;
      const std::array<double, 2> c{0.5 * (p[0] + q[0]) - nudge * size * ey / len,
                                    0.5 * (p[1] + q[1]) + nudge * size * ex / len};
      if (f.contains(c[0], c[1])) return c;
    }
  return std::nullopt;
}

// Distance of a parameter point from the loops of a face (polyline approximation).
double distance_to_loops(const TrimmedSurface& f, double u, double v) {
  double d = 1e300;
  const auto seg = [&](const std::vector<std::array<double, 2>>& poly) {
    for (std::size_t i = 0, j = poly.size() - 1; i < poly.size(); j = i++) {
      const double ax = poly[j][0], ay = poly[j][1], bx = poly[i][0], by = poly[i][1];
      const double ex = bx - ax, ey = by - ay, l2 = ex * ex + ey * ey;
      const double t = l2 > 0.0 ? std::clamp(((u - ax) * ex + (v - ay) * ey) / l2, 0.0, 1.0) : 0.0;
      d = std::min(d, std::hypot(u - ax - t * ex, v - ay - t * ey));
    }
  };
  seg(f.outer().polyline(16));
  for (const auto& h : f.holes()) seg(h.polyline(16));
  return d;
}

// Parameter intervals of an intersection curve on which both p-curves lie inside their faces.
std::vector<std::pair<double, double>> shared_intervals(const IntersectionCurve& c, const TrimmedSurface& fa,
                                                        const TrimmedSurface& fb) {
  const auto inside = [&](double t) {
    const Vec3 w1 = c.pcurve1.evaluate(t), w2 = c.pcurve2.evaluate(t);
    return fa.contains(w1.x, w1.y) && fb.contains(w2.x, w2.y);
  };
  const double t0 = c.curve.domain_start(), t1 = c.curve.domain_end();
  const int n = 256;
  std::vector<std::pair<double, double>> out;
  bool prev = inside(t0);
  double start = t0, tp = t0;
  const auto edge = [&](double lo, double hi, bool lo_in) {
    for (int it = 0; it < 50; ++it) {
      const double m = 0.5 * (lo + hi);
      (inside(m) == lo_in ? lo : hi) = m;
    }
    return 0.5 * (lo + hi);
  };
  for (int i = 1; i <= n; ++i) {
    const double t = t0 + (t1 - t0) * i / n;
    const bool cur = inside(t);
    if (cur != prev) {
      const double e = edge(tp, t, prev);
      if (prev) out.emplace_back(start, e);
      else start = e;
      prev = cur;
    }
    tp = t;
  }
  if (prev) out.emplace_back(start, t1);
  const double tiny = 1e-9 * (t1 - t0);
  out.erase(std::remove_if(out.begin(), out.end(), [&](const auto& iv) { return iv.second - iv.first <= tiny; }),
            out.end());
  return out;
}

}  // namespace

bool point_in_shell(const Shell& shell, const Vec3& p, double tolerance) {
  const Box3 b = shell_box(shell);
  const double L = 10.0 * (distance(b.lo, b.hi) + distance(p, 0.5 * (b.lo + b.hi)) + 1.0);
  const Vec3 dirs[] = {normalized(Vec3{0.5377, 0.8013, 0.2622}), normalized(Vec3{-0.3112, 0.4419, -0.8414}),
                       normalized(Vec3{0.7071, -0.6018, 0.3712}), normalized(Vec3{-0.8123, -0.2311, 0.5354}),
                       normalized(Vec3{0.1234, -0.9512, -0.2832})};
  int votes_in = 0, votes = 0;
  for (const auto& d : dirs) {
    const NurbsCurve ray = make_line(p, p + L * d);
    int count = 0;
    bool valid = true;
    for (const auto& f : shell.faces()) {
      for (const auto& h : intersect(ray, f.surface(), {tolerance * 1e-3, 1e-10, 1e-9})) {
        if (!f.contains(h.u, h.v)) continue;
        const auto du = f.surface().domain_u(), dv = f.surface().domain_v();
        const double psize = std::max(du[1] - du[0], dv[1] - dv[0]);
        Vec3 n;
        try {
          n = f.surface().normal(h.u, h.v);
        } catch (const GeometryError&) {
          valid = false;
          break;
        }
        if (distance_to_loops(f, h.u, h.v) < 1e-6 * psize || std::abs(dot(n, d)) < 1e-3) {
          valid = false;  // boundary hit or grazing ray: ambiguous
          break;
        }
        ++count;
      }
      if (!valid) break;
    }
    if (!valid) continue;
    ++votes;
    if (count % 2 == 1) ++votes_in;
    if (votes == 3) break;
  }
  if (votes == 0) throw GeometryError("point_in_shell: no unambiguous ray found");
  return 2 * votes_in > votes;
}

Shell boolean(const Shell& a_in, const Shell& b_in, BooleanOp op, double tolerance) {
  if (!a_in.is_closed() || !b_in.is_closed()) throw GeometryError("boolean: both shells must be closed");
  const Shell A = a_in.oriented(), B = b_in.oriented();
  const auto& FA = A.faces();
  const auto& FB = B.faces();
  std::vector<std::vector<NurbsCurve>> cuts_a(FA.size()), cuts_b(FB.size());
  SurfaceIntersectionOptions opt;
  opt.tol.length = 1e-3 * tolerance;
  for (std::size_t i = 0; i < FA.size(); ++i)
    for (std::size_t j = 0; j < FB.size(); ++j) {
      if (!face_box(FA[i]).overlaps(face_box(FB[j]), tolerance)) continue;
      const auto res = intersect(FA[i].surface(), FB[j].surface(), opt);
      for (const auto& c : res.curves)
        for (const auto& [t0, t1] : shared_intervals(c, FA[i], FB[j])) {
          cuts_a[i].push_back(c.pcurve1.sub_curve(t0, t1));
          cuts_b[j].push_back(c.pcurve2.sub_curve(t0, t1));
        }
    }
  const auto pieces = [&](const std::vector<TrimmedSurface>& faces, const std::vector<std::vector<NurbsCurve>>& cuts) {
    std::vector<TrimmedSurface> out;
    for (std::size_t i = 0; i < faces.size(); ++i) {
      const auto parts = split(faces[i], cuts[i]);
      out.insert(out.end(), parts.begin(), parts.end());
    }
    return out;
  };
  const auto inside_of = [&](const TrimmedSurface& piece, const Shell& other) {
    const auto q = interior_point(piece);
    if (!q) throw GeometryError("boolean: cannot find an interior point of a face piece");
    return point_in_shell(other, piece.surface().evaluate((*q)[0], (*q)[1]), tolerance);
  };
  std::vector<TrimmedSurface> result;
  for (const auto& p : pieces(FA, cuts_a)) {
    const bool in = inside_of(p, B);
    if ((op == BooleanOp::Intersection) == in) result.push_back(p);
  }
  for (const auto& p : pieces(FB, cuts_b)) {
    const bool in = inside_of(p, A);
    if (op == BooleanOp::Union && !in) result.push_back(p);
    if (op == BooleanOp::Intersection && in) result.push_back(p);
    if (op == BooleanOp::Difference && in) result.push_back(p.flipped());
  }
  // the pieces keep the orientation of their (outward oriented) shells; difference pieces of B
  // are flipped – re-orienting here would turn inner voids outwards
  return Shell::sew(std::move(result), tolerance);
}

}  // namespace nf
