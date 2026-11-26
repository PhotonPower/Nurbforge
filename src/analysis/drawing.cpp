#include "nurbforge/analysis/drawing.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <limits>
#include <optional>
#include <string>
#include <utility>

#include "nurbforge/analysis/closest_point.hpp"
#include "nurbforge/analysis/intersect.hpp"
#include "nurbforge/core/error.hpp"
#include "nurbforge/geometry/primitives.hpp"
#include "nurbforge/geometry/transform.hpp"

namespace nf {

namespace {

struct Box {
  Vec3 lo{1e300, 1e300, 1e300}, hi{-1e300, -1e300, -1e300};
  void add(const Vec3& p) {
    lo = {std::min(lo.x, p.x), std::min(lo.y, p.y), std::min(lo.z, p.z)};
    hi = {std::max(hi.x, p.x), std::max(hi.y, p.y), std::max(hi.z, p.z)};
  }
  [[nodiscard]] bool contains(const Vec3& p, double eps) const {
    return p.x >= lo.x - eps && p.x <= hi.x + eps && p.y >= lo.y - eps && p.y <= hi.y + eps &&
           p.z >= lo.z - eps && p.z <= hi.z + eps;
  }
  [[nodiscard]] bool overlaps(const Box& o, double eps) const {
    return lo.x <= o.hi.x + eps && o.lo.x <= hi.x + eps && lo.y <= o.hi.y + eps &&
           o.lo.y <= hi.y + eps && lo.z <= o.hi.z + eps && o.lo.z <= hi.z + eps;
  }
  // slab test of the segment P + t w, t in [0, L], against the box grown by eps
  [[nodiscard]] bool hit_by(const Vec3& P, const Vec3& w, double L, double eps) const {
    double t0 = 0.0, t1 = L;
    const double p[3] = {P.x, P.y, P.z}, d[3] = {w.x, w.y, w.z};
    const double a[3] = {lo.x - eps, lo.y - eps, lo.z - eps};
    const double b[3] = {hi.x + eps, hi.y + eps, hi.z + eps};
    for (int k = 0; k < 3; ++k) {
      if (std::abs(d[k]) < 1e-300) {
        if (p[k] < a[k] || p[k] > b[k]) return false;
        continue;
      }
      double ta = (a[k] - p[k]) / d[k], tb = (b[k] - p[k]) / d[k];
      if (ta > tb) std::swap(ta, tb);
      t0 = std::max(t0, ta);
      t1 = std::min(t1, tb);
      if (t0 > t1) return false;
    }
    return true;
  }
};

Box box_of(const NurbsCurve& c) {
  Box b;
  for (const auto& p : c.control_points()) b.add(p);
  return b;
}

// The faces of the occluding shells with their bounding boxes.
class Occluders {
 public:
  explicit Occluders(const std::vector<Shell>& shells) {
    for (const auto& sh : shells) {
      const std::size_t first = faces_.size();
      for (const auto& f : sh.faces()) {
        Box b;
        for (const auto& p : f.surface().control_points()) {
          b.add(p);
          all_.add(p);
        }
        faces_.push_back({&f, b, {}});
      }
      for (const auto& e : sh.edges())
        for (const auto& use : e.uses) faces_[first + use.face].edges.push_back(&e.curve);
    }
    size_ = faces_.empty() ? 1.0 : distance(all_.lo, all_.hi);
  }

  [[nodiscard]] double size() const { return size_; }

  // Distance from P to beyond the occluders in the direction w.
  [[nodiscard]] double reach(const Vec3& P, const Vec3& w) const {
    if (faces_.empty()) return 0.0;
    double L = 0.0;
    for (int k = 0; k < 8; ++k) {
      const Vec3 c{k & 1 ? all_.hi.x : all_.lo.x, k & 2 ? all_.hi.y : all_.lo.y,
                   k & 4 ? all_.hi.z : all_.lo.z};
      L = std::max(L, dot(c - P, w));
    }
    return L + 1e-3 * size_;
  }

  // True if a face is crossed by the ray from P to the viewer: the point lies more than tol off
  // the tangent plane of the face where the ray meets it.
  [[nodiscard]] bool hides(const Vec3& P, const View& view, double tol) const {
    const Vec3 w = view.to_viewer(P);
    const double L = view.eye ? distance(*view.eye, P) : reach(P, w);
    if (!(L > tol)) return false;
    const NurbsCurve ray = make_line(P, P + L * w);
    const Tolerance itol{std::min(1e-9, 1e-3 * tol), tol::kParametric, tol::kAngular};
    for (const auto& face : faces_) {
      if (!face.box.hit_by(P, w, L, tol)) continue;
      const TrimmedSurface* f = face.face;
      for (const auto& h : intersect(ray, f->surface(), itol)) {
        const double t = dot(h.point - P, w);
        if (t <= tol || !face.has(h, tol)) continue;
        try {
          if (t * std::abs(dot(f->surface().normal(h.u, h.v), w)) > tol) return true;
        } catch (const GeometryError&) {
          // no normal at the hit (degenerate point): not counted
        }
      }
    }
    return false;
  }

  // True if the line through P along the unit vector d meets a face.
  [[nodiscard]] bool covers(const Vec3& P, const Vec3& d) const {
    const double fwd = reach(P, d), back = reach(P, -d);
    if (!(fwd + back > 0.0)) return false;
    const Vec3 A = P - back * d;
    const double L = fwd + back;
    const NurbsCurve ray = make_line(A, A + L * d);
    for (const auto& face : faces_) {
      if (!face.box.hit_by(A, d, L, 1e-9 * size_)) continue;
      for (const auto& h : intersect(ray, face.face->surface()))
        if (face.face->contains(h.u, h.v)) return true;
    }
    return false;
  }

 private:
  struct Face {
    const TrimmedSurface* face;
    Box box;
    std::vector<const NurbsCurve*> edges;  // 3D boundary curves (edges of the shell)

    // The hit lies in the trimmed region or, within eps, on its boundary. A ray through an
    // edge of the shell (a hidden edge of a cube exactly behind a visible one) is on the
    // boundary of both faces, where contains() may give either answer for each of them.
    [[nodiscard]] bool has(const CurveSurfaceIntersection& h, double eps) const {
      if (face->contains(h.u, h.v)) return true;
      for (const NurbsCurve* e : edges)
        if (box_of(*e).contains(h.point, eps) && closest_point(*e, h.point).distance <= eps)
          return true;
      return false;
    }
  };
  std::vector<Face> faces_;
  Box all_;
  double size_{1.0};
};

struct Interval {
  double t0{0.0}, t1{0.0};
  bool state{false};
};

// Maximal intervals of [a, b] with a constant value of the predicate, from n + 1 samples; the
// changes are refined by bisection.
std::vector<Interval> intervals(const std::function<bool(double)>& pred, double a, double b,
                                std::size_t n) {
  std::vector<Interval> out;
  double ta = a, start = a;
  bool sa = pred(a);
  for (std::size_t k = 1; k <= n; ++k) {
    const double tb = k == n ? b : a + (b - a) * static_cast<double>(k) / static_cast<double>(n);
    const bool sb = pred(tb);
    if (sb != sa) {
      double lo = ta, hi = tb;
      for (int it = 0; it < 60 && hi - lo > 1e-10 * (b - a); ++it) {
        const double m = 0.5 * (lo + hi);
        (pred(m) == sa ? lo : hi) = m;
      }
      const double x = 0.5 * (lo + hi);
      out.push_back({start, x, sa});
      start = x;
    }
    ta = tb;
    sa = sb;
  }
  out.push_back({start, b, sa});
  return out;
}

NurbsCurve part(const NurbsCurve& c, double t0, double t1) {
  const double a = c.domain_start(), b = c.domain_end(), eps = 1e-9 * (b - a);
  const bool at_a = t0 <= a + eps, at_b = t1 >= b - eps;
  return at_a && at_b ? c : c.sub_curve(at_a ? a : t0, at_b ? b : t1);
}

std::size_t span_count(const NurbsCurve& c) {
  return std::max<std::size_t>(1, unique_knots(c.knots()).first.size() - 1);
}

// Drawing plane of a view: through the origin, normal to the view direction (parallel) or to
// the line from the origin to the eye (perspective).
struct Frame {
  Vec3 origin, s, t, n;
  std::optional<Vec3> eye;
  double eye_distance{0.0};  // of the eye from the plane

  // Projection of a curve onto the plane, in plane coordinates (s, t, 0). Parallel: an affine
  // map of the control points. Perspective: a projective map, exact too: the control point P
  // with weight w goes to d ((P - E) . s, (P - E) . t) / z with the weight w z, where z is its
  // depth in front of the eye E and d the distance of the eye from the plane.
  [[nodiscard]] NurbsCurve project(const NurbsCurve& c) const {
    if (!eye)
      return c.transformed(Transform(
          {s.x, s.y, s.z, -dot(s, origin), t.x, t.y, t.z, -dot(t, origin), 0.0, 0.0, 0.0, 0.0}));
    std::vector<Vec3> P;
    std::vector<double> W;
    P.reserve(c.control_points().size());
    for (std::size_t i = 0; i < c.control_points().size(); ++i) {
      const Vec3 r = c.control_points()[i] - *eye;
      const double z = -dot(r, n);
      if (!(z > 1e-9 * eye_distance))
        throw GeometryError("perspective view: the model reaches behind the eye point");
      P.push_back({eye_distance * dot(r, s) / z, eye_distance * dot(r, t) / z, 0.0});
      W.push_back(c.weights()[i] * z);
    }
    return {c.degree(), c.knots(), std::move(P), std::move(W)};
  }
  [[nodiscard]] Vec3 to_space(const Vec3& q) const { return origin + q.x * s + q.y * t; }
};

Frame frame_for(const View& view, const DrawingOptions& opt, const char* who) {
  if (!(opt.tolerance > 0.0))
    throw GeometryError(std::string(who) + ": tolerance must be positive");
  Frame f;
  f.origin = opt.origin;
  if (view.eye) {
    f.eye = view.eye;
    f.eye_distance = distance(*view.eye, opt.origin);
    if (!(f.eye_distance > 0.0))
      throw GeometryError(std::string(who) +
                          ": the eye point coincides with the origin of the drawing plane");
    f.n = (*view.eye - opt.origin) / f.eye_distance;
  } else {
    const double len = norm(view.direction);
    if (!(len > 0.0) || !std::isfinite(len))
      throw GeometryError(std::string(who) + ": the direction must be a non-zero vector");
    f.n = -1.0 * view.direction / len;
  }
  Vec3 s = opt.s_axis - dot(opt.s_axis, f.n) * f.n;
  if (!(norm(s) > 1e-9 * std::max(1.0, norm(opt.s_axis)))) {
    if (norm(opt.s_axis) > 0.0)
      throw GeometryError(std::string(who) + ": s_axis is parallel to the view direction");
    const Vec3 helper = std::abs(f.n.z) < 0.9 ? Vec3{0, 0, 1} : Vec3{0, 1, 0};
    s = cross(helper, f.n);  // horizontal for views that are not vertical
  }
  f.s = normalized(s);
  f.t = cross(f.n, f.s);
  return f;
}

double join_tolerance(const DrawingOptions& opt) {
  return opt.join_tolerance > 0.0 ? opt.join_tolerance : 10.0 * opt.tolerance;
}

// A curve of the model with its source.
struct Item {
  NurbsCurve c3;
  DrawingCurveKind kind{DrawingCurveKind::Edge};
  std::size_t shell{0}, index{0};
};

// Parameters where the projection of the curve folds back onto itself (its tangent reverses:
// the 3D tangent passes through the view direction, e.g. a circle seen edge-on).
std::vector<double> fold_parameters(const NurbsCurve& c2) {
  const double a = c2.domain_start(), b = c2.domain_end();
  const std::size_t n = std::max<std::size_t>(64, 8 * span_count(c2));
  std::vector<double> out;
  double ta = a;
  Vec3 da = c2.derivatives(a, 1)[1];
  double scale = 0.0;  // largest speed so far: at a fold the speed vanishes
  for (std::size_t k = 1; k <= n; ++k) {
    const double tb = k == n ? b : a + (b - a) * static_cast<double>(k) / static_cast<double>(n);
    const Vec3 db = c2.derivatives(tb, 1)[1];
    const double na = norm(da), nb = norm(db);
    scale = std::max({scale, na, nb});
    if (!(nb > 1e-9 * scale)) continue;  // a sample on a fold: compare across it
    if (!(na > 1e-9 * scale)) {
      ta = tb;
      da = db;
      continue;
    }
    if (dot(da, db) < 0.0 && norm(cross(da, db)) < 0.5 * na * nb) {
      double lo = ta, hi = tb;
      for (int it = 0; it < 60 && hi - lo > 1e-12 * (b - a); ++it) {
        const double m = 0.5 * (lo + hi);
        (dot(c2.derivatives(m, 1)[1], da) > 0.0 ? lo : hi) = m;
      }
      out.push_back(0.5 * (lo + hi));
    }
    ta = tb;
    da = db;
  }
  return out;
}

// True if the two faces at the edge have the same tangent plane along it (the seam of a closed
// face, a fillet running into a face): such an edge is no line of a drawing.
bool is_smooth_edge(const Shell& sh, const ShellEdge& e) {
  if (e.uses.size() != 2) return false;
  const double a = e.curve.domain_start(), b = e.curve.domain_end();
  for (const double fr : {0.2, 0.5, 0.8}) {
    const Vec3 P = e.curve.evaluate(a + fr * (b - a));
    Vec3 n[2];
    for (int k = 0; k < 2; ++k) {
      const Coedge& use = e.uses[static_cast<std::size_t>(k)];
      const TrimmedSurface& f = sh.faces()[use.face];
      const TrimLoop& loop = use.loop == 0 ? f.outer() : f.holes()[use.loop - 1];
      const NurbsCurve& pc = loop.curves[use.index];
      const double tau = use.reversed ? 1.0 - fr : fr;
      const Vec3 uv = pc.evaluate(pc.domain_start() + tau * (pc.domain_end() - pc.domain_start()));
      try {
        const SurfacePoint sp = closest_point_local(f.surface(), P, uv.x, uv.y);
        n[k] = f.surface().normal(sp.u, sp.v);
      } catch (const GeometryError&) {
        return false;
      }
    }
    if (norm(cross(n[0], n[1])) > 1e-6) return false;
  }
  return true;
}

// Edges and silhouettes of the shells, split where their projections fold.
std::vector<Item> model_curves(const std::vector<Shell>& shells, const View& view,
                               const Frame& frame, const DrawingOptions& opt) {
  std::vector<Item> raw;
  SilhouetteOptions so;
  so.tolerance = opt.tolerance;
  so.resolution = opt.resolution;
  for (std::size_t k = 0; k < shells.size(); ++k) {
    const auto& edges = shells[k].edges();
    for (std::size_t e = 0; e < edges.size(); ++e)
      if (opt.tangent_edges || !is_smooth_edge(shells[k], edges[e]))
        raw.push_back({edges[e].curve, DrawingCurveKind::Edge, k, e});
    for (auto& sc : silhouette(shells[k], view, so))
      raw.push_back({std::move(sc.curve), DrawingCurveKind::Silhouette, k, sc.face});
  }
  std::vector<Item> out;
  for (auto& it : raw) {
    const auto folds = fold_parameters(frame.project(it.c3));
    if (folds.empty()) {
      out.push_back(std::move(it));
      continue;
    }
    std::vector<double> cuts{it.c3.domain_start()};
    const double len = it.c3.domain_end() - it.c3.domain_start();
    for (const double t : folds)
      if (t - cuts.back() > 1e-9 * len && it.c3.domain_end() - t > 1e-9 * len) cuts.push_back(t);
    cuts.push_back(it.c3.domain_end());
    for (std::size_t i = 0; i + 1 < cuts.size(); ++i)
      out.push_back({part(it.c3, cuts[i], cuts[i + 1]), it.kind, it.shell, it.index});
  }
  return out;
}

double extent_of(const NurbsCurve& c) {
  const double a = c.domain_start(), b = c.domain_end();
  const Vec3 first = c.evaluate(a);
  double e = 0.0;
  for (int i = 1; i <= 16; ++i)
    e = std::max(e, distance(c.evaluate(a + (b - a) * i / 16.0), first));
  return e;
}

// Lines already drawn, for the removal of duplicates.
class Drawn {
 public:
  explicit Drawn(double jt) : jt_(jt) {}

  void add(const NurbsCurve& c2) {
    Entry e{c2, box_of(c2), {}, {}};
    const double a = c2.domain_start(), b = c2.domain_end();
    const std::size_t n = 16 * span_count(c2);
    for (std::size_t i = 0; i <= n; ++i) {
      e.params.push_back(i == n ? b
                                : a + (b - a) * static_cast<double>(i) / static_cast<double>(n));
      e.points.push_back(c2.evaluate(e.params.back()));
    }
    curves_.push_back(std::move(e));
  }

  // The parts of c2 that do not lie on a line already drawn (parameter intervals). Stretches on
  // another line shorter than a few tolerances (crossings, common end points) are kept.
  [[nodiscard]] std::vector<std::array<double, 2>> uncovered(const NurbsCurve& c2) const {
    const double a = c2.domain_start(), b = c2.domain_end();
    const Box box = box_of(c2);
    std::vector<const Entry*> near;
    for (const auto& e : curves_)
      if (e.box.overlaps(box, jt_)) near.push_back(&e);
    if (near.empty()) return {{a, b}};
    const auto covered = [&](double t) {
      const auto d = c2.derivatives(t, 1);
      const double speed = norm(d[1]);
      for (const Entry* e : near) {
        if (!e->box.contains(d[0], jt_)) continue;
        const auto [u, dist] = e->nearest(d[0]);
        if (dist > jt_) continue;
        // beyond the end of the other line (its continuation): not on it
        const double range = e->params.back() - e->params.front();
        if ((u - e->params.front() < 1e-9 * range || e->params.back() - u < 1e-9 * range) &&
            dist > 1e-3 * jt_)
          continue;
        const Vec3 tk = e->curve.derivatives(u, 1)[1];
        const double nk = norm(tk);
        if (!(speed > 0.0) || !(nk > 0.0) || norm(cross(d[1], tk)) <= 0.05 * speed * nk)
          return true;
      }
      return false;
    };
    auto parts = intervals(covered, a, b, std::max<std::size_t>(32, 8 * span_count(c2)));
    for (auto& p : parts)
      if (p.state && distance(c2.evaluate(p.t0), c2.evaluate(p.t1)) < 20.0 * jt_ &&
          !(p.t0 <= a && p.t1 >= b))
        p.state = false;
    std::vector<std::array<double, 2>> out;
    for (const auto& p : parts) {
      if (p.state) continue;
      if (!out.empty() && out.back()[1] == p.t0)
        out.back()[1] = p.t1;
      else
        out.push_back({p.t0, p.t1});
    }
    std::erase_if(out, [&](const auto& iv) { return !(iv[1] - iv[0] > 1e-7 * (b - a)); });
    return out;
  }

 private:
  struct Entry {
    NurbsCurve curve;
    Box box;
    std::vector<double> params;  // samples for the nearest-point search
    std::vector<Vec3> points;

    // Parameter and distance of the point of the curve nearest to p: the nearest sample,
    // refined by a golden-section search between its neighbours. No derivatives are used, as
    // a projected curve may stop where it folds (zero speed), where Newton's method fails.
    [[nodiscard]] std::pair<double, double> nearest(const Vec3& p) const {
      std::size_t best = 0;
      for (std::size_t i = 1; i < points.size(); ++i)
        if (distance(points[i], p) < distance(points[best], p)) best = i;
      double lo = params[best == 0 ? 0 : best - 1];
      double hi = params[std::min(best + 1, params.size() - 1)];
      const auto f = [&](double t) { return distance(curve.evaluate(t), p); };
      constexpr double g = 0.6180339887498949;
      double x1 = hi - g * (hi - lo), x2 = lo + g * (hi - lo), f1 = f(x1), f2 = f(x2);
      for (int it = 0; it < 60 && hi - lo > 1e-13 * (params.back() - params.front()); ++it) {
        if (f1 < f2) {
          hi = x2, x2 = x1, f2 = f1;
          x1 = hi - g * (hi - lo), f1 = f(x1);
        } else {
          lo = x1, x1 = x2, f1 = f2;
          x2 = lo + g * (hi - lo), f2 = f(x2);
        }
      }
      const double t = 0.5 * (lo + hi);
      return {t, f(t)};
    }
  };
  double jt_;
  std::vector<Entry> curves_;
};

}  // namespace

namespace {

// Visibility of the curves. `outlines` are the edges and silhouettes of the occluders: the
// visibility of a curve can change only where its projection crosses the projection of one of
// them (or where the curve passes through a face). Each curve is cut at these crossings, found
// by intersecting the projected curves; within the intervals in between the visibility is
// sampled and changes are refined by bisection, which finds the passages through faces and
// guards against a missed silhouette. `same_index`: outlines[k] is curves[k] itself.
std::vector<CurvePiece> visibility_of(const std::vector<NurbsCurve>& curves, const View& view,
                                      const Occluders& occ, const std::vector<NurbsCurve>& outlines,
                                      bool same_index, const VisibilityOptions& opt) {
  // projections for the crossings, in a plane in front of the eye through the middle of
  // everything; without one (a curve behind the eye) only the sampling is used
  std::vector<NurbsCurve> curves2, outlines2;
  std::vector<Box> boxes2;
  if (opt.events && !outlines.empty()) {
    Box all;
    for (const auto& c : curves)
      for (const auto& p : c.control_points()) all.add(p);
    for (const auto& c : outlines)
      for (const auto& p : c.control_points()) all.add(p);
    DrawingOptions dopt;
    dopt.tolerance = opt.tolerance;
    dopt.origin = 0.5 * (all.lo + all.hi);
    try {
      const Frame frame = frame_for(view, dopt, "visibility");
      for (const auto& c : curves) curves2.push_back(frame.project(c));
      for (const auto& c : outlines) {
        outlines2.push_back(frame.project(c));
        boxes2.push_back(box_of(outlines2.back()));
      }
    } catch (const GeometryError&) {
      curves2.clear();
      outlines2.clear();
    }
  }
  const Tolerance itol{std::max(opt.tolerance, 1e-12), tol::kParametric, tol::kAngular};

  std::vector<CurvePiece> out;
  for (std::size_t k = 0; k < curves.size(); ++k) {
    const NurbsCurve& c = curves[k];
    const double a = c.domain_start(), b = c.domain_end();
    const auto visible = [&](double t) { return !occ.hides(c.evaluate(t), view, opt.tolerance); };
    // cuts: the crossings with the projected outlines
    std::vector<double> cuts{a};
    if (!curves2.empty()) {
      std::vector<double> x;
      const Box box = box_of(curves2[k]);
      for (std::size_t j = 0; j < outlines2.size(); ++j) {
        if ((same_index && j == k) || !box.overlaps(boxes2[j], opt.tolerance)) continue;
        for (const auto& h : intersect(curves2[k], outlines2[j], itol)) x.push_back(h.s);
      }
      std::sort(x.begin(), x.end());
      for (const double t : x)
        if (t - cuts.back() > 1e-9 * (b - a) && b - t > 1e-9 * (b - a)) cuts.push_back(t);
    }
    cuts.push_back(b);
    // sampling between the cuts, at interior points (the state at a crossing is ambiguous)
    const std::size_t n = std::max<std::size_t>({opt.samples, 4 * span_count(c), 2});
    std::vector<Interval> raw;
    for (std::size_t i = 0; i + 1 < cuts.size(); ++i) {
      const double t0 = cuts[i], t1 = cuts[i + 1];
      const auto m = std::max<std::size_t>(
          1, static_cast<std::size_t>(std::ceil(static_cast<double>(n) * (t1 - t0) / (b - a))));
      const double h = (t1 - t0) / static_cast<double>(m);
      const double first = t0 + 0.5 * h, last = t1 - 0.5 * h;
      std::vector<Interval> part_ivs = m == 1 ? std::vector<Interval>{{first, last, visible(first)}}
                                              : intervals(visible, first, last, m - 1);
      part_ivs.front().t0 = t0;
      part_ivs.back().t1 = t1;
      raw.insert(raw.end(), part_ivs.begin(), part_ivs.end());
    }
    // stretches shorter than a few tolerances (the end point of a hidden edge on a visible
    // one) are no pieces of their own
    std::vector<Interval> ivs;
    for (const auto& iv : raw) {
      const bool tiny = distance(c.evaluate(iv.t0), c.evaluate(iv.t1)) <= 10.0 * opt.tolerance &&
                        iv.t1 - iv.t0 < 0.5 * (b - a) / static_cast<double>(n);
      if (!ivs.empty() && (tiny || ivs.back().state == iv.state))
        ivs.back().t1 = iv.t1;
      else
        ivs.push_back(iv);
    }
    if (ivs.size() > 1 && distance(c.evaluate(ivs[0].t0), c.evaluate(ivs[0].t1)) <=
                              10.0 * opt.tolerance) {  // tiny first piece
      ivs[1].t0 = ivs[0].t0;
      ivs.erase(ivs.begin());
      if (ivs.size() > 1 && ivs[0].state == ivs[1].state) {
        ivs[1].t0 = ivs[0].t0;
        ivs.erase(ivs.begin());
      }
    }
    for (const auto& iv : ivs) out.push_back({part(c, iv.t0, iv.t1), iv.state, k, iv.t0, iv.t1});
  }
  return out;
}

}  // namespace

std::vector<CurvePiece> visibility(const std::vector<NurbsCurve>& curves, const View& view,
                                   const std::vector<Shell>& occluders,
                                   const VisibilityOptions& opt) {
  if (!(opt.tolerance > 0.0)) throw GeometryError("visibility: tolerance must be positive");
  const Occluders occ(occluders);
  std::vector<NurbsCurve> outlines;
  if (opt.events) {
    SilhouetteOptions so;
    so.tolerance = opt.tolerance;
    so.resolution = opt.resolution;
    for (const auto& sh : occluders) {
      for (const auto& e : sh.edges()) outlines.push_back(e.curve);
      for (auto& sc : silhouette(sh, view, so)) outlines.push_back(std::move(sc.curve));
    }
  }
  return visibility_of(curves, view, occ, outlines, false, opt);
}

Drawing make2d(const std::vector<Shell>& shells, const Vec3& direction, const DrawingOptions& opt) {
  return make2d(shells, View::parallel(direction), opt);
}

Drawing make2d(const std::vector<Shell>& shells, const View& view, const DrawingOptions& opt) {
  const Frame frame = frame_for(view, opt, "make2d");
  const double jt = join_tolerance(opt);
  const std::vector<Item> items = model_curves(shells, view, frame, opt);
  std::vector<NurbsCurve> curves;
  curves.reserve(items.size());
  for (const auto& it : items) curves.push_back(it.c3);
  VisibilityOptions vo;
  vo.tolerance = opt.tolerance;
  vo.samples = opt.samples;
  vo.resolution = opt.resolution;
  // the model curves are the outlines of the occluders themselves
  const std::vector<CurvePiece> pieces =
      visibility_of(curves, view, Occluders(shells), curves, true, vo);

  Drawing d{frame.origin, frame.s, frame.t, frame.n, {}, {}};
  Drawn drawn(jt);
  for (const bool want_visible : {true, false}) {
    if (!want_visible && !opt.hidden) break;
    for (const auto& p : pieces) {
      if (p.visible != want_visible) continue;
      const NurbsCurve c2 = frame.project(p.curve);
      if (extent_of(c2) <= jt) continue;  // seen end-on
      const Item& src = items[p.source];
      for (const auto& iv : drawn.uncovered(c2)) {
        NurbsCurve piece2 = part(c2, iv[0], iv[1]);
        if (extent_of(piece2) <= jt) continue;
        drawn.add(piece2);
        (want_visible ? d.visible : d.hidden)
            .push_back({std::move(piece2), part(p.curve, iv[0], iv[1]), want_visible, src.kind,
                        src.shell, src.index});
      }
    }
  }
  return d;
}

Outline outline(const std::vector<Shell>& shells, const Vec3& direction,
                const DrawingOptions& opt) {
  return outline(shells, View::parallel(direction), opt);
}

Outline outline(const std::vector<Shell>& shells, const View& view, const DrawingOptions& opt) {
  const Frame frame = frame_for(view, opt, "outline");
  const double jt = join_tolerance(opt);
  const Occluders occ(shells);

  // the distinct projected lines
  std::vector<NurbsCurve> lines;
  {
    Drawn drawn(jt);
    for (const auto& it : model_curves(shells, view, frame, opt)) {
      const NurbsCurve c2 = frame.project(it.c3);
      if (extent_of(c2) <= jt) continue;
      for (const auto& iv : drawn.uncovered(c2)) {
        NurbsCurve piece = part(c2, iv[0], iv[1]);
        if (extent_of(piece) <= jt) continue;
        drawn.add(piece);
        lines.push_back(std::move(piece));
      }
    }
  }
  // split at their mutual intersections
  std::vector<std::vector<double>> cuts(lines.size());
  {
    std::vector<Box> boxes;
    for (const auto& c : lines) boxes.push_back(box_of(c));
    const Tolerance itol{jt, tol::kParametric, tol::kAngular};
    for (std::size_t i = 0; i < lines.size(); ++i)
      for (std::size_t j = i + 1; j < lines.size(); ++j) {
        if (!boxes[i].overlaps(boxes[j], jt)) continue;
        for (const auto& h : intersect(lines[i], lines[j], itol)) {
          cuts[i].push_back(h.s);
          cuts[j].push_back(h.t);
        }
      }
  }
  // the pieces with the shells on exactly one side, oriented with that side on the left
  std::vector<NurbsCurve> pieces;
  const double probe = 20.0 * jt;
  for (std::size_t i = 0; i < lines.size(); ++i) {
    const NurbsCurve& c = lines[i];
    std::vector<double> t{c.domain_start()};
    std::sort(cuts[i].begin(), cuts[i].end());
    for (const double x : cuts[i])
      if (distance(c.evaluate(x), c.evaluate(t.back())) > jt &&
          distance(c.evaluate(x), c.evaluate(c.domain_end())) > jt)
        t.push_back(x);
    t.push_back(c.domain_end());
    for (std::size_t k = 0; k + 1 < t.size(); ++k) {
      const double tm = 0.5 * (t[k] + t[k + 1]);
      const auto dm = c.derivatives(tm, 1);
      const double speed = norm(dm[1]);
      if (!(speed > 0.0)) continue;
      const Vec3 left = Vec3{-dm[1].y, dm[1].x, 0.0} / speed;
      // the view ray through a point of the drawing plane
      const auto covered = [&](const Vec3& q) {
        const Vec3 p = frame.to_space(q);
        return occ.covers(p, frame.eye ? normalized(p - *frame.eye) : frame.n);
      };
      const bool on_left = covered(dm[0] + probe * left);
      const bool on_right = covered(dm[0] - probe * left);
      if (on_left == on_right) continue;
      NurbsCurve piece = part(c, t[k], t[k + 1]);
      pieces.push_back(on_left ? std::move(piece) : piece.reversed());
    }
  }
  // chain into closed loops
  Outline out{frame.origin, frame.s, frame.t, frame.n, {}};
  const double ct = 4.0 * jt;
  const auto start_of = [](const NurbsCurve& c) { return c.evaluate(c.domain_start()); };
  const auto end_of = [](const NurbsCurve& c) { return c.evaluate(c.domain_end()); };
  std::vector<char> used(pieces.size(), 0);
  for (std::size_t i = 0; i < pieces.size(); ++i) {
    if (used[i]) continue;
    used[i] = 1;
    std::vector<NurbsCurve> loop{pieces[i]};
    const Vec3 first = start_of(pieces[i]);
    Vec3 at = end_of(pieces[i]);
    while (distance(at, first) > ct || (loop.size() == 1 && extent_of(loop[0]) <= ct)) {
      std::size_t best = pieces.size();
      double dbest = ct;
      for (std::size_t j = 0; j < pieces.size(); ++j) {
        if (used[j]) continue;
        const double dj = distance(start_of(pieces[j]), at);
        if (dj <= dbest) {
          dbest = dj;
          best = j;
        }
      }
      if (best == pieces.size())
        throw GeometryError("outline: the contour does not close at (s, t) = (" +
                            std::to_string(at.x) + ", " + std::to_string(at.y) + ")");
      used[best] = 1;
      at = end_of(pieces[best]);
      loop.push_back(pieces[best]);
    }
    out.loops.push_back(std::move(loop));
  }
  return out;
}

}  // namespace nf
