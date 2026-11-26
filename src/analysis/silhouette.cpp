#include "nurbforge/analysis/silhouette.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <limits>
#include <map>
#include <numbers>
#include <optional>
#include <string>
#include <utility>

#include "clip_detail.hpp"
#include "nurbforge/analysis/closest_point.hpp"
#include "nurbforge/core/error.hpp"
#include "nurbforge/fitting/fit.hpp"
#include "nurbforge/geometry/primitives.hpp"
#include "nurbforge/modeling/construct.hpp"
#include "nurbforge/nurbs/join.hpp"

namespace nf {

View View::parallel(const Vec3& direction) {
  View v;
  v.direction = direction;
  return v;
}

View View::perspective(const Vec3& eye) {
  View v;
  v.eye = eye;
  return v;
}

Vec3 View::to_viewer(const Vec3& p) const {
  const Vec3 w = eye ? *eye - p : -direction;
  const double n = norm(w);
  if (!(n > 0.0) || !std::isfinite(n))
    throw GeometryError(eye ? "view: a point coincides with the eye point"
                            : "view: the direction must be a non-zero vector");
  return w / n;
}

namespace {

std::vector<double> spans_in(const KnotVector& U, double a, double b, std::size_t sub) {
  std::vector<double> x{a};
  for (const double k : unique_knots(U).first)
    if (k > a + tol::kParametric && k < b - tol::kParametric) x.push_back(k);
  x.push_back(b);
  std::vector<double> out{a};
  for (std::size_t i = 0; i + 1 < x.size(); ++i)
    for (std::size_t q = 1; q <= sub; ++q)
      out.push_back(q == sub ? x[i + 1]
                             : x[i] + (x[i + 1] - x[i]) * static_cast<double>(q) /
                                          static_cast<double>(sub));
  return out;
}

// Interior knots across which the surface is only C0 (multiplicity >= degree).
std::vector<double> c0_knots(const KnotVector& U, std::size_t p, double a, double b) {
  std::vector<double> out;
  const auto [k, m] = unique_knots(U);
  for (std::size_t i = 0; i < k.size(); ++i)
    if (k[i] > a + tol::kParametric && k[i] < b - tol::kParametric &&
        static_cast<std::size_t>(m[i]) >= p)
      out.push_back(k[i]);
  return out;
}

// g = N . w on a smooth piece of the surface.
class Facing {
 public:
  Facing(const NurbsSurface& s, const View& view) : s_(s), view_(view) {}
  [[nodiscard]] double operator()(double u, double v) const {
    return dot(s_.normal(u, v), view_.to_viewer(s_.evaluate(u, v)));
  }

 private:
  const NurbsSurface& s_;
  const View& view_;
};

// Zero set of a function on a grid (marching squares). A value counts as positive unless it is
// clearly negative, so a zero set running exactly along grid lines or through grid nodes (a
// silhouette on a knot line of the surface) is traced along them without ambiguity.
struct Polyline {
  std::vector<std::array<double, 2>> points;
  bool closed{false};
};

std::vector<Polyline> zero_set(const std::function<double(double, double)>& f,
                               const std::vector<double>& xs, const std::vector<double>& ys,
                               double& max_abs) {
  const std::size_t nx = xs.size() - 1, ny = ys.size() - 1;
  const auto positive = [](double v) { return v > -1e-12; };
  std::vector<char> pos((nx + 1) * (ny + 1));
  max_abs = 0.0;
  bool any_pos = false, any_neg = false;
  for (std::size_t i = 0; i <= nx; ++i)
    for (std::size_t j = 0; j <= ny; ++j) {
      const double v = f(xs[i], ys[j]);
      max_abs = std::max(max_abs, std::abs(v));
      pos[i * (ny + 1) + j] = positive(v) ? 1 : 0;
      (positive(v) ? any_pos : any_neg) = true;
    }
  if (!any_pos || !any_neg) return {};
  const auto P = [&](std::size_t i, std::size_t j) { return pos[i * (ny + 1) + j] != 0; };
  // edges: horizontal (i, j)-(i+1, j) with id j * nx + i, vertical (i, j)-(i, j+1) after them
  const std::size_t nh = nx * (ny + 1);
  const auto hid = [&](std::size_t i, std::size_t j) { return j * nx + i; };
  const auto vid = [&](std::size_t i, std::size_t j) { return nh + j * (nx + 1) + i; };
  const std::size_t nedges = nh + (nx + 1) * ny;
  constexpr std::size_t kNone = static_cast<std::size_t>(-1);
  std::vector<std::array<std::size_t, 2>> link(nedges, {kNone, kNone});
  const auto connect = [&](std::size_t a, std::size_t b) {
    (link[a][0] == kNone ? link[a][0] : link[a][1]) = b;
    (link[b][0] == kNone ? link[b][0] : link[b][1]) = a;
  };
  for (std::size_t i = 0; i < nx; ++i)
    for (std::size_t j = 0; j < ny; ++j) {
      const bool s00 = P(i, j), s10 = P(i + 1, j), s01 = P(i, j + 1), s11 = P(i + 1, j + 1);
      const std::size_t bottom = hid(i, j), top = hid(i, j + 1), left = vid(i, j),
                        right = vid(i + 1, j);
      std::vector<std::size_t> e;
      if (s00 != s10) e.push_back(bottom);
      if (s10 != s11) e.push_back(right);
      if (s01 != s11) e.push_back(top);
      if (s00 != s01) e.push_back(left);
      if (e.size() == 2) {
        connect(e[0], e[1]);
      } else if (e.size() == 4) {  // saddle: the centre decides which corners are cut off
        const bool c = positive(f(0.5 * (xs[i] + xs[i + 1]), 0.5 * (ys[j] + ys[j + 1])));
        if (c == s00) {
          connect(bottom, right);
          connect(left, top);
        } else {
          connect(bottom, left);
          connect(top, right);
        }
      }
    }
  // the crossing on an edge: where f changes from negative to positive (bisection)
  const auto crossing = [&](std::size_t id) {
    double ax, ay, bx, by;
    if (id < nh) {
      const std::size_t j = id / nx, i = id % nx;
      ax = xs[i], bx = xs[i + 1], ay = by = ys[j];
    } else {
      const std::size_t r = id - nh, j = r / (nx + 1), i = r % (nx + 1);
      ax = bx = xs[i], ay = ys[j], by = ys[j + 1];
    }
    const bool pa = positive(f(ax, ay));
    double lo = 0.0, hi = 1.0;
    for (int it = 0; it < 48; ++it) {
      const double m = 0.5 * (lo + hi);
      (positive(f(ax + m * (bx - ax), ay + m * (by - ay))) == pa ? lo : hi) = m;
    }
    const double m = pa ? lo : hi;  // on the non-negative side
    return std::array<double, 2>{ax + m * (bx - ax), ay + m * (by - ay)};
  };
  std::vector<char> used(nedges, 0);
  std::vector<Polyline> out;
  const auto walk = [&](std::size_t start, bool closed) {
    Polyline pl;
    pl.closed = closed;
    std::size_t prev = kNone, cur = start;
    while (cur != kNone && !used[cur]) {
      used[cur] = 1;
      const auto q = crossing(cur);
      if (pl.points.empty() ||
          std::hypot(q[0] - pl.points.back()[0], q[1] - pl.points.back()[1]) > 1e-12)
        pl.points.push_back(q);
      const std::size_t next = link[cur][0] != prev ? link[cur][0] : link[cur][1];
      prev = cur;
      cur = next;
    }
    if (closed && pl.points.size() > 2) pl.points.push_back(pl.points.front());
    if (pl.points.size() >= 2) out.push_back(std::move(pl));
  };
  for (std::size_t e = 0; e < nedges; ++e)  // open chains start on the boundary
    if (!used[e] && link[e][0] != kNone && link[e][1] == kNone) walk(e, false);
  for (std::size_t e = 0; e < nedges; ++e)
    if (!used[e] && link[e][0] != kNone) walk(e, true);
  return out;
}

// Cubic interpolation of a family of curves sharing one parameter: `exact(t)` gives one point
// per curve, `error(curves, t, exact(t))` the deviation at t. Parameters are added at the
// midpoints of failing intervals.
std::vector<NurbsCurve> fit_family(
    std::vector<double> T, const std::function<std::vector<Vec3>(double)>& exact,
    const std::function<double(const std::vector<NurbsCurve>&, double, const std::vector<Vec3>&)>&
        error,
    double tol) {
  std::map<double, std::vector<Vec3>> cache;
  const auto at = [&](double t) -> const std::vector<Vec3>& {
    auto it = cache.find(t);
    if (it == cache.end()) it = cache.emplace(t, exact(t)).first;
    return it->second;
  };
  const std::size_t ncurves = at(T.front()).size();
  for (int round = 0; round < 16; ++round) {
    std::vector<std::vector<Vec3>> pts(ncurves);
    for (const double t : T)
      for (std::size_t k = 0; k < ncurves; ++k) pts[k].push_back(at(t)[k]);
    const std::size_t deg = std::min<std::size_t>(3, T.size() - 1);
    std::vector<NurbsCurve> curves;
    for (std::size_t k = 0; k < ncurves; ++k) curves.push_back(interpolate_curve(pts[k], deg, T));
    std::vector<double> add;
    for (std::size_t i = 0; i + 1 < T.size(); ++i) {
      const double tm = 0.5 * (T[i] + T[i + 1]);
      if (error(curves, tm, at(tm)) > tol) add.push_back(tm);
    }
    if (add.empty()) return curves;
    if (T.size() + add.size() > 8193) break;
    T.insert(T.end(), add.begin(), add.end());
    std::sort(T.begin(), T.end());
  }
  throw GeometryError("silhouette: tolerance " + std::to_string(tol) + " not reached");
}

// Parameters where f changes sign along [a, b] (n samples, bisection).
std::vector<double> sign_changes(const std::function<double(double)>& f, double a, double b,
                                 std::size_t n) {
  std::vector<double> out;
  double ta = a, fa = f(a);
  for (std::size_t k = 1; k <= n; ++k) {
    const double tb = k == n ? b : a + (b - a) * static_cast<double>(k) / static_cast<double>(n);
    const double fb = f(tb);
    if (fb == 0.0) continue;  // wait for a definite sign
    if (fa != 0.0 && (fa > 0.0) != (fb > 0.0)) {
      double lo = ta, hi = tb, flo = fa;
      for (int it = 0; it < 60 && hi - lo > 1e-13 * (b - a); ++it) {
        const double mid = 0.5 * (lo + hi), fm = f(mid);
        if (fm == 0.0) {
          lo = hi = mid;
        } else if ((fm > 0.0) == (flo > 0.0)) {
          lo = mid;
          flo = fm;
        } else {
          hi = mid;
        }
      }
      out.push_back(0.5 * (lo + hi));
    }
    ta = tb;
    fa = fb;
  }
  return out;
}

void find_cusps(SilhouetteCurve& c, const NurbsSurface& s, const View& view) {
  const double a = c.curve.domain_start(), b = c.curve.domain_end();
  const auto f = [&](double t) {
    const auto d = c.curve.derivatives(t, 1);
    const double speed = norm(d[1]);
    if (!(speed > 0.0)) return 0.0;
    const Vec3 q = c.pcurve.evaluate(t);
    const double val = dot(cross(d[1] / speed, view.to_viewer(d[0])), s.normal(q.x, q.y));
    return std::abs(val) < 1e-9 ? 0.0 : val;
  };
  const std::size_t spans = unique_knots(c.curve.knots()).first.size();
  c.cusps = sign_changes(f, a, b, std::max<std::size_t>(64, 4 * spans));
}

// Smooth silhouettes of a piece of the surface on which g is continuous.
void smooth_silhouettes(const NurbsSurface& s, const View& view, const SilhouetteOptions& opt,
                        std::vector<SilhouetteCurve>& out) {
  const auto du = s.domain_u(), dv = s.domain_v();
  const double lu = du[1] - du[0], lv = dv[1] - dv[0];
  const Facing g(s, view);

  // Sides of the parameter domain that collapse to a point (poles) and on which g vanishes (a
  // view direction normal to the axis of a sphere): such a side would be part of the zero set
  // and a branch of every silhouette through the pole, and next to it the normal is too
  // inaccurate for the sign of g. The grid stops a little before such sides; its curves are
  // stretched to the whole domain, so those ending there reach the pole (the strip left out is
  // so narrow in 3D that the position of the curve within it does not matter).
  double size = 0.0;
  for (const auto& q : s.control_points())
    size = std::max(size, distance(q, s.control_points().front()));
  const auto collapsed = [&](bool u_side, bool at_end) {
    const std::size_t nu = s.num_u(), nv = s.num_v();
    const Vec3 ref = s.control_point(u_side && at_end ? nu - 1 : 0, !u_side && at_end ? nv - 1 : 0);
    double extent = 0.0;
    for (std::size_t k = 0; k < (u_side ? nv : nu); ++k)
      extent = std::max(extent, distance(ref, u_side ? s.control_point(at_end ? nu - 1 : 0, k)
                                                     : s.control_point(k, at_end ? nv - 1 : 0)));
    if (extent > 1e-12 * std::max(size, 1e-300)) return false;
    // only where g vanishes on that side; otherwise no silhouette reaches the pole
    for (int k = 0; k <= 8; ++k) {
      const double w = k / 8.0;
      const double val = u_side ? g(at_end ? du[1] : du[0], dv[0] + w * lv)
                                : g(du[0] + w * lu, at_end ? dv[1] : dv[0]);
      if (std::abs(val) > 1e-9) return false;
    }
    return true;
  };
  const bool pole[4] = {collapsed(true, false), collapsed(true, true), collapsed(false, false),
                        collapsed(false, true)};
  constexpr double eps = 1e-5;
  // the domain of the grid in normalised parameters
  const double x0 = pole[0] ? eps : 0.0, x1 = pole[1] ? 1.0 - eps : 1.0;
  const double y0 = pole[2] ? eps : 0.0, y1 = pole[3] ? 1.0 - eps : 1.0;
  // g in normalised parameters, divided by the distances to those sides: its zero set no
  // longer contains them, and its gradient near a pole points across the silhouette
  const auto scaled = [&](double x, double y) {
    double w = 1.0;
    if (pole[0]) w *= std::max(x, 0.5 * eps);
    if (pole[1]) w *= std::max(1.0 - x, 0.5 * eps);
    if (pole[2]) w *= std::max(y, 0.5 * eps);
    if (pole[3]) w *= std::max(1.0 - y, 0.5 * eps);
    return g(du[0] + x * lu, dv[0] + y * lv) / w;
  };
  // grid: the knot spans of the surface, each divided into opt.resolution cells
  const auto grid = [&](const KnotVector& U, double a, double b, double lo, double len) {
    std::vector<double> x;
    for (const double k : spans_in(U, a, b, std::max<std::size_t>(opt.resolution, 2)))
      x.push_back((k - lo) / len);
    return x;
  };
  const auto xs = grid(s.knots_u(), du[0] + x0 * lu, du[0] + x1 * lu, du[0], lu);
  const auto ys = grid(s.knots_v(), dv[0] + y0 * lv, dv[0] + y1 * lv, dv[0], lv);
  double max_abs = 0.0;
  const auto lines = zero_set(scaled, xs, ys, max_abs);
  if (max_abs < 1e-9) return;  // seen edge-on everywhere

  // exact silhouette point next to a point (normalised parameters x, y): Newton steps along the
  // gradient; points on the boundary of the grid stay on it
  const auto polish = [&](double x, double y) {
    if (x < x0 || x > x1 || y < y0 || y > y1) return Vec3{du[0] + x * lu, dv[0] + y * lv, 0.0};
    const bool fix_x = x <= x0 + 1e-9 || x >= x1 - 1e-9, fix_y = y <= y0 + 1e-9 || y >= y1 - 1e-9;
    const auto G = [&](double a, double b) {
      return scaled(std::clamp(a, x0, x1), std::clamp(b, y0, y1));
    };
    constexpr double hstep = 1e-6;
    for (int it = 0; it < 8 && !(fix_x && fix_y); ++it) {
      const double val = G(x, y);
      if (std::abs(val) < 1e-14) break;
      const double xa = std::max(x0, x - hstep), xb = std::min(x1, x + hstep);
      const double ya = std::max(y0, y - hstep), yb = std::min(y1, y + hstep);
      const double gx = fix_x ? 0.0 : (G(xb, y) - G(xa, y)) / (xb - xa);
      const double gy = fix_y ? 0.0 : (G(x, yb) - G(x, ya)) / (yb - ya);
      const double g2 = gx * gx + gy * gy;
      // differences of g carry round-off of about 1e-9: a smaller gradient is noise
      if (!(g2 > 1e-12)) break;
      const double sx = -val * gx / g2, sy = -val * gy / g2;
      const double len = std::hypot(sx, sy);
      // The traced point is much closer to the zero set than a cell of the grid. A long step
      // means that g does not vary in the free directions (a point held on the boundary, with
      // the zero set running along it): the point stays where it is.
      if (len > 0.02) break;
      x = std::clamp(x + sx, x0, x1);
      y = std::clamp(y + sy, y0, y1);
      if (len < 1e-15) break;
    }
    return Vec3{du[0] + x * lu, dv[0] + y * lv, 0.0};
  };
  std::vector<double> knots_x, knots_y;  // interior knot lines in normalised parameters
  for (const double k : unique_knots(s.knots_u()).first)
    if (k > du[0] + (x0 + 1e-9) * lu && k < du[0] + (x1 - 1e-9) * lu)
      knots_x.push_back((k - du[0]) / lu);
  for (const double k : unique_knots(s.knots_v()).first)
    if (k > dv[0] + (y0 + 1e-9) * lv && k < dv[0] + (y1 - 1e-9) * lv)
      knots_y.push_back((k - dv[0]) / lv);

  for (const auto& line : lines) {
    // smooth guide curve through the crossings of the grid lines (normalised parameters)
    std::vector<Vec3> pts;
    std::vector<double> par;
    for (const auto& q : line.points) {
      const Vec3 p{q[0], q[1], 0.0};
      if (!pts.empty())
        par.push_back(par.back() + distance(p, pts.back()));
      else
        par.push_back(0.0);
      pts.push_back(p);
    }
    const NurbsCurve trace = interpolate_curve(pts, std::min<std::size_t>(3, pts.size() - 1), par);
    const double a = trace.domain_start(), b = trace.domain_end();
    const auto exact = [&](double t) {
      const Vec3 q = trace.evaluate(t);
      // the grid's domain stretched to the whole parameter domain
      const Vec3 uv = polish(std::clamp((q.x - x0) / (x1 - x0), 0.0, 1.0),
                             std::clamp((q.y - y0) / (y1 - y0), 0.0, 1.0));
      return std::vector<Vec3>{s.evaluate(uv.x, uv.y), uv};
    };
    const auto error = [&](const std::vector<NurbsCurve>& c, double t, const std::vector<Vec3>& e) {
      const Vec3 q = c[1].evaluate(t);
      return std::max(
          distance(c[0].evaluate(t), e[0]),
          distance(s.evaluate(std::clamp(q.x, du[0], du[1]), std::clamp(q.y, dv[0], dv[1])), e[0]));
    };
    // The curve is fitted piecewise between the parameters where it crosses a knot line of the
    // surface (it is smooth in between only); the pieces are joined.
    const std::size_t n = std::max<std::size_t>(64, 4 * unique_knots(trace.knots()).first.size());
    std::vector<double> cuts{a, b};
    for (const double k : knots_x)
      for (const double t :
           sign_changes([&](double w) { return trace.evaluate(w).x - k; }, a, b, n))
        cuts.push_back(t);
    for (const double k : knots_y)
      for (const double t :
           sign_changes([&](double w) { return trace.evaluate(w).y - k; }, a, b, n))
        cuts.push_back(t);
    std::sort(cuts.begin(), cuts.end());
    // cuts next to each other or to the ends are merged
    {
      std::vector<double> merged{a};
      for (const double c : cuts)
        if (c - merged.back() > 1e-6 * (b - a) && b - c > 1e-6 * (b - a)) merged.push_back(c);
      merged.push_back(b);
      cuts = std::move(merged);
    }

    std::optional<SilhouetteCurve> chain;
    for (std::size_t i = 0; i + 1 < cuts.size(); ++i) {
      const double t0 = cuts[i], t1 = cuts[i + 1];
      std::vector<double> T;
      for (int k = 0; k <= 4; ++k) T.push_back(k == 4 ? t1 : t0 + (t1 - t0) * k / 4.0);
      auto fit = fit_family(T, exact, error, opt.tolerance);
      if (!chain) {
        chain = SilhouetteCurve{std::move(fit[0]), std::move(fit[1]), false, false, 0, {}};
      } else {
        chain->curve = join(chain->curve, fit[0], opt.tolerance);
        chain->pcurve = join(chain->pcurve, fit[1], 1e-6 * std::max(lu, lv));
      }
    }
    if (!chain) continue;
    // a curve that does not move in 3D (along a collapsed side) is no silhouette
    double extent = 0.0;
    const Vec3 first = chain->curve.evaluate(a);
    for (int k = 1; k <= 16; ++k)
      extent = std::max(extent, distance(chain->curve.evaluate(a + (b - a) * k / 16.0), first));
    if (extent <= opt.tolerance) continue;
    chain->closed = distance(first, chain->curve.evaluate(b)) <= opt.tolerance;
    find_cusps(*chain, s, view);
    out.push_back(std::move(*chain));
  }
}

// Parts of the knot line (u = x, or v = x if !along_u... ) across which g changes sign.
void crease_silhouettes(const NurbsSurface& s, const NurbsSurface& lo, const NurbsSurface& hi,
                        bool u_line, double x, const View& view,
                        std::vector<SilhouetteCurve>& out) {
  const auto dw = u_line ? s.domain_v() : s.domain_u();
  const auto side = [&](const NurbsSurface& part, double w) {
    const double val = u_line ? Facing(part, view)(x, w) : Facing(part, view)(w, x);
    return std::abs(val) < 1e-9 ? 0.0 : val;
  };
  // +1 where the two sides face opposite ways, -1 where they agree
  const auto f = [&](double w) {
    const double a = side(lo, w), b = side(hi, w);
    if (a == 0.0 || b == 0.0) return 0.0;
    return (a > 0.0) != (b > 0.0) ? 1.0 : -1.0;
  };
  std::vector<double> cuts{dw[0]};
  for (const double w : sign_changes(f, dw[0], dw[1], 64)) cuts.push_back(w);
  cuts.push_back(dw[1]);
  for (std::size_t i = 0; i + 1 < cuts.size(); ++i) {
    const double a = cuts[i], b = cuts[i + 1];
    if (!(b - a > 1e-9 * (dw[1] - dw[0])) || !(f(0.5 * (a + b)) > 0.0)) continue;
    const NurbsCurve iso = u_line ? s.iso_curve_u(x) : s.iso_curve_v(x);
    const Vec3 p0 = u_line ? Vec3{x, a, 0} : Vec3{a, x, 0},
               p1 = u_line ? Vec3{x, b, 0} : Vec3{b, x, 0};
    SilhouetteCurve sc{
        iso.sub_curve(a, b), reparametrize(make_line(p0, p1), a, b), false, true, 0, {}};
    out.push_back(std::move(sc));
  }
}

// True if g jumps across the knot line (one-sided values from the two parts).
bool is_crease(const NurbsSurface& lo, const NurbsSurface& hi, bool u_line, double x,
               const std::array<double, 2>& dw, const View& view) {
  const Facing a(lo, view), b(hi, view);
  for (int i = 0; i <= 16; ++i) {
    const double w = dw[0] + (dw[1] - dw[0]) * (i + 0.5) / 17.0;
    const double d = u_line ? a(x, w) - b(x, w) : a(w, x) - b(w, x);
    if (std::abs(d) > 1e-8) return true;
  }
  return false;
}

}  // namespace

namespace {

// Removes curves lying on an earlier one (a silhouette on the seam of a closed surface is found
// on both sides of the parameter domain).
void remove_duplicates(std::vector<SilhouetteCurve>& curves, double tol) {
  const auto covered = [&](const NurbsCurve& c, const NurbsCurve& by) {
    const double a = c.domain_start(), b = c.domain_end();
    for (int i = 0; i <= 8; ++i)
      if (closest_point(by, c.evaluate(a + (b - a) * i / 8.0)).distance > 4.0 * tol) return false;
    return true;
  };
  std::vector<SilhouetteCurve> out;
  for (auto& c : curves) {
    bool dup = false;
    for (const auto& o : out)
      if (covered(c.curve, o.curve) && covered(o.curve, c.curve)) {
        dup = true;
        break;
      }
    if (!dup) out.push_back(std::move(c));
  }
  curves = std::move(out);
}

}  // namespace

std::vector<SilhouetteCurve> silhouette(const NurbsSurface& surface, const View& view,
                                        const SilhouetteOptions& opt) {
  if (!(opt.tolerance > 0.0)) throw GeometryError("silhouette: tolerance must be positive");
  (void)view.to_viewer(surface.control_points().front() + Vec3{1e-3, 2e-3, 3e-3});
  const NurbsSurface s = surface.clamped();
  const auto du = s.domain_u(), dv = s.domain_v();
  std::vector<SilhouetteCurve> out;

  // knot lines across which the normal jumps split the surface into smooth pieces
  std::vector<double> cu{du[0]}, cv{dv[0]};
  for (const double x : c0_knots(s.knots_u(), s.degree_u(), du[0], du[1])) {
    const auto [lo, hi] = s.split_u(x);
    if (!is_crease(lo, hi, true, x, dv, view)) continue;
    cu.push_back(x);
    crease_silhouettes(s, lo, hi, true, x, view, out);
  }
  for (const double x : c0_knots(s.knots_v(), s.degree_v(), dv[0], dv[1])) {
    const auto [lo, hi] = s.split_v(x);
    if (!is_crease(lo, hi, false, x, du, view)) continue;
    cv.push_back(x);
    crease_silhouettes(s, lo, hi, false, x, view, out);
  }
  cu.push_back(du[1]);
  cv.push_back(dv[1]);
  if (cu.size() == 2 && cv.size() == 2) {
    smooth_silhouettes(s, view, opt, out);
  } else {
    for (std::size_t i = 0; i + 1 < cu.size(); ++i)
      for (std::size_t j = 0; j + 1 < cv.size(); ++j)
        smooth_silhouettes(s.sub_surface(cu[i], cu[i + 1], cv[j], cv[j + 1]), view, opt, out);
  }
  remove_duplicates(out, opt.tolerance);
  return out;
}

std::vector<SilhouetteCurve> silhouette(const TrimmedSurface& face, const View& view,
                                        const SilhouetteOptions& opt) {
  std::vector<SilhouetteCurve> all = silhouette(face.surface(), view, opt);
  if (face.is_untrimmed()) return all;
  std::vector<SilhouetteCurve> out;
  for (auto& c : all) {
    const double a = c.pcurve.domain_start(), b = c.pcurve.domain_end();
    for (const auto& iv : detail::inside_intervals(face, c.pcurve)) {
      if (!(iv[1] - iv[0] > 1e-9 * (b - a))) continue;
      if (iv[0] <= a && iv[1] >= b) {
        out.push_back(c);
        continue;
      }
      SilhouetteCurve part{c.curve.sub_curve(iv[0], iv[1]),
                           c.pcurve.sub_curve(iv[0], iv[1]),
                           false,
                           c.crease,
                           c.face,
                           {}};
      for (const double t : c.cusps)
        if (t > iv[0] && t < iv[1]) part.cusps.push_back(t);
      out.push_back(std::move(part));
    }
  }
  return out;
}

std::vector<SilhouetteCurve> silhouette(const Shell& shell, const View& view,
                                        const SilhouetteOptions& opt) {
  std::vector<SilhouetteCurve> out;
  for (std::size_t f = 0; f < shell.faces().size(); ++f)
    for (auto& c : silhouette(shell.faces()[f], view, opt)) {
      c.face = f;
      out.push_back(std::move(c));
    }
  return out;
}

}  // namespace nf
