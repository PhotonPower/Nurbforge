#include "nurbforge/analysis/shell_section.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <limits>
#include <optional>
#include <string>

#include "nurbforge/core/error.hpp"
#include "nurbforge/geometry/transform.hpp"

namespace nf {

namespace {

// ------------------------------------------------------------------ helpers

double model_size(const Shell& sh) {
  Vec3 lo{1e300, 1e300, 1e300}, hi{-1e300, -1e300, -1e300};
  for (const auto& f : sh.faces())
    for (const auto& p : f.surface().control_points()) {
      lo = {std::min(lo.x, p.x), std::min(lo.y, p.y), std::min(lo.z, p.z)};
      hi = {std::max(hi.x, p.x), std::max(hi.y, p.y), std::max(hi.z, p.z)};
    }
  return std::max(distance(lo, hi), 1e-300);
}

/// Degree-1 curve (u, v, 0) from a to b on the parameter interval [t0, t1].
NurbsCurve line_on(const Vec3& a, const Vec3& b, double t0, double t1) {
  return NurbsCurve(1, {t0, t0, t1, t1}, {a, b});
}

/// If the height n . S(u, v) of the surface depends on v only (every row j of control points has
/// one height c_j and the weights factor as w_ij = a_i b_j): the height as a function of v, i.e.
/// the 1D NURBS with control points (c_j, 0, 0) and weights b_j. Transposed for u.
std::optional<NurbsCurve> height_curve(const NurbsSurface& s, const Vec3& n, bool along_v,
                                       double tol) {
  const std::size_t nu = s.num_u(), nv = s.num_v();
  const std::size_t rows = along_v ? nv : nu, cols = along_v ? nu : nv;
  const auto P = [&](std::size_t r, std::size_t c) {
    return along_v ? s.control_point(c, r) : s.control_point(r, c);
  };
  const auto W = [&](std::size_t r, std::size_t c) {
    return along_v ? s.weight(c, r) : s.weight(r, c);
  };
  std::vector<Vec3> h(rows);
  std::vector<double> b(rows);
  for (std::size_t r = 0; r < rows; ++r) {
    const double c0 = dot(n, P(r, 0));
    for (std::size_t c = 1; c < cols; ++c)
      if (std::abs(dot(n, P(r, c)) - c0) > tol) return std::nullopt;
    h[r] = {c0, 0, 0};
    b[r] = W(r, 0);
  }
  // product form of the weights: w(r, c) / w(r, 0) independent of r
  for (std::size_t c = 1; c < cols; ++c) {
    const double q = W(0, c) / W(0, 0);
    for (std::size_t r = 1; r < rows; ++r)
      if (std::abs(W(r, c) / W(r, 0) - q) > 1e-12 * std::abs(q)) return std::nullopt;
  }
  return NurbsCurve(along_v ? s.degree_v() : s.degree_u(), along_v ? s.knots_v() : s.knots_u(),
                    std::move(h), std::move(b));
}

/// Roots of h(t) = c on the domain of the height curve (sign changes of a dense sampling, refined
/// by bisection to the resolution of double).
std::vector<double> height_roots(const NurbsCurve& h, double c) {
  const double t0 = h.domain_start(), t1 = h.domain_end();
  const std::size_t spans = std::max<std::size_t>(1, unique_knots(h.knots()).first.size() - 1);
  const std::size_t n = 64 * spans;
  const auto f = [&](double t) { return h.evaluate(t).x - c; };
  std::vector<double> roots;
  double ta = t0, fa = f(t0);
  if (fa == 0.0) roots.push_back(t0);
  for (std::size_t k = 1; k <= n; ++k) {
    const double tb =
        k == n ? t1 : t0 + (t1 - t0) * static_cast<double>(k) / static_cast<double>(n);
    const double fb = f(tb);
    if (fb == 0.0) {
      roots.push_back(tb);
    } else if ((fa < 0.0) != (fb < 0.0) && fa != 0.0) {
      double lo = ta, hi = tb, flo = fa;
      for (int it = 0; it < 200 && hi - lo > 4 * std::numeric_limits<double>::epsilon() *
                                                 std::max(1.0, std::abs(lo));
           ++it) {
        const double mid = 0.5 * (lo + hi), fm = f(mid);
        if ((fm < 0.0) == (flo < 0.0)) {
          lo = mid;
          flo = fm;
        } else {
          hi = mid;
        }
      }
      roots.push_back(0.5 * (lo + hi));
    }
    ta = tb;
    fa = fb;
  }
  return roots;
}

/// Raw section curves of one face surface (3D curve and p-curve with the same parameter).
struct Piece {
  NurbsCurve c3, pc;
  bool exact;
};

std::vector<Piece> surface_section(const NurbsSurface& s, const Plane& plane, const Vec3& n,
                                   double size, const SectionOptions& opt) {
  const double c = dot(n, plane.point), tol = 1e-12 * size;
  const auto hv = height_curve(s, n, true, tol), hu = height_curve(s, n, false, tol);
  std::vector<Piece> out;
  if (hv && hu) {  // constant height: a face parallel to the plane
    if (std::abs(hv->control_points()[0].x - c) <= 1e-9 * size)
      throw GeometryError(
          "section: the plane contains a face (choose a height between critical heights)");
    return out;
  }
  const auto du = s.domain_u(), dv = s.domain_v();
  if (hv) {
    for (double v : height_roots(*hv, c))
      out.push_back({s.iso_curve_v(v), line_on({du[0], v, 0}, {du[1], v, 0}, du[0], du[1]), true});
    return out;
  }
  if (hu) {
    for (double u : height_roots(*hu, c))
      out.push_back({s.iso_curve_u(u), line_on({u, dv[0], 0}, {u, dv[1], 0}, dv[0], dv[1]), true});
    return out;
  }
  SurfaceIntersectionOptions so;
  so.fit_tolerance = opt.fit_tolerance > 0.0 ? opt.fit_tolerance : 1e-9 * size;
  for (auto& ic : planar_section(s, plane, so).curves)
    out.push_back({std::move(ic.curve), std::move(ic.pcurve1), false});
  return out;
}

/// Parameter intervals of the p-curve inside the trimmed region (bisection on contains()).
std::vector<std::array<double, 2>> inside_intervals(const TrimmedSurface& f, const NurbsCurve& pc) {
  const double t0 = pc.domain_start(), t1 = pc.domain_end();
  if (f.is_untrimmed()) return {{t0, t1}};
  const auto in = [&](double t) {
    const Vec3 p = pc.evaluate(t);
    return f.contains(p.x, p.y);
  };
  const std::size_t spans = std::max<std::size_t>(1, unique_knots(pc.knots()).first.size() - 1);
  const std::size_t n = std::max<std::size_t>(128, 32 * spans);
  const auto refine = [&](double a, double b, bool ina) {  // a has state ina, b the other state
    for (int it = 0; it < 200 && b - a > 4 * std::numeric_limits<double>::epsilon() *
                                             std::max(1.0, std::abs(a));
         ++it) {
      const double m = 0.5 * (a + b);
      (in(m) == ina ? a : b) = m;
    }
    return 0.5 * (a + b);
  };
  std::vector<std::array<double, 2>> out;
  double ta = t0;
  bool ina = in(t0);
  double start = ina ? t0 : 0.0;
  for (std::size_t k = 1; k <= n; ++k) {
    const double tb =
        k == n ? t1 : t0 + (t1 - t0) * static_cast<double>(k) / static_cast<double>(n);
    const bool inb = in(tb);
    if (inb != ina) {
      const double x = refine(ta, tb, ina);
      if (ina)
        out.push_back({start, x});
      else
        start = x;
    }
    ta = tb;
    ina = inb;
  }
  if (ina) out.push_back({start, t1});
  return out;
}

NurbsCurve to_plane(const NurbsCurve& c3, const Vec3& o, const Vec3& e1, const Vec3& e2) {
  return c3.transformed(
      Transform({e1.x, e1.y, e1.z, -dot(e1, o), e2.x, e2.y, e2.z, -dot(e2, o), 0, 0, 0, 0}));
}

/// Signed area 1/2 sum int (s t' - t s') of a loop in plane coordinates (8-point Gauss per span).
double loop_area(const std::vector<SectionCurve>& curves) {
  static const double gx[8] = {-0.9602898564975363, -0.7966664774136267, -0.5255324099163290,
                               -0.1834346424956498, 0.1834346424956498,  0.5255324099163290,
                               0.7966664774136267,  0.9602898564975363};
  static const double gw[8] = {0.1012285362903763, 0.2223810344533745, 0.3137066458778873,
                               0.3626837833783620, 0.3626837833783620, 0.3137066458778873,
                               0.2223810344533745, 0.1012285362903763};
  double a = 0.0;
  for (const auto& sc : curves) {
    const auto k = unique_knots(sc.curve.knots()).first;
    for (std::size_t i = 0; i + 1 < k.size(); ++i) {
      if (k[i + 1] <= k[i] || k[i + 1] <= sc.curve.domain_start() || k[i] >= sc.curve.domain_end())
        continue;
      // 4 sub-intervals per span (rational pieces, e.g. quarter circles)
      for (int q = 0; q < 4; ++q) {
        const double a0 = k[i] + (k[i + 1] - k[i]) * q / 4.0,
                     a1 = k[i] + (k[i + 1] - k[i]) * (q + 1) / 4.0;
        for (int g = 0; g < 8; ++g) {
          const double t = 0.5 * (a0 + a1) + 0.5 * (a1 - a0) * gx[g];
          const auto d = sc.curve.derivatives(t, 1);
          a += 0.25 * (a1 - a0) * gw[g] * (d[0].x * d[1].y - d[0].y * d[1].x);
        }
      }
    }
  }
  return a;
}

}  // namespace

Vec3 SectionCurve::normal(double t) const {
  const Vec3 p = pcurve.evaluate(t);
  return normal_sign * surface.normal(p.x, p.y);
}

PlanarSection section(const Shell& input, const Plane& plane, const SectionOptions& opt) {
  const Shell sh = input.is_consistently_oriented() ? input : input.oriented();
  const double size = model_size(sh);
  const Vec3 n = normalized(plane.normal);
  PlanarSection out;
  out.origin = plane.point;
  out.normal = n;
  {
    Vec3 s = opt.s_axis - dot(opt.s_axis, n) * n;
    if (norm(s) < 1e-12) {
      const Vec3 helper = std::abs(n.x) < 0.9 ? Vec3{1, 0, 0} : Vec3{0, 1, 0};
      s = cross(n, helper);
    }
    out.s_axis = normalized(s);
    out.t_axis = cross(n, out.s_axis);
  }
  // outward orientation: a closed, consistently oriented shell has normals pointing outwards if its
  // volume is positive
  const double sign = sh.is_closed() && sh.volume(1e-6) < 0.0 ? -1.0 : 1.0;
  // pieces of all faces, clipped to the trimmed regions
  std::vector<SectionCurve> pieces;
  const double min_len = 1e-3 * opt.join_tolerance;
  for (std::size_t f = 0; f < sh.faces().size(); ++f) {
    const TrimmedSurface& face = sh.faces()[f];
    for (const Piece& p : surface_section(face.surface(), plane, n, size, opt))
      for (const auto& iv : inside_intervals(face, p.pc)) {
        const bool whole = iv[0] == p.pc.domain_start() && iv[1] == p.pc.domain_end();
        NurbsCurve c3 = whole ? p.c3 : p.c3.sub_curve(iv[0], iv[1]);
        if (distance(c3.evaluate(c3.domain_start()), c3.evaluate(c3.domain_end())) < min_len &&
            norm(c3.control_points().front() - c3.control_points()[c3.num_control_points() / 2]) <
                min_len)
          continue;  // degenerate piece (pole, touching trim)
        NurbsCurve pc = whole ? p.pc : p.pc.sub_curve(iv[0], iv[1]);
        // curve (plane coordinates) is set after chaining
        pieces.push_back(SectionCurve{c3, c3, std::move(pc), f, p.exact, face.surface(), sign});
      }
  }
  // chain the pieces into closed loops
  const auto start = [](const SectionCurve& c) {
    return c.curve3d.evaluate(c.curve3d.domain_start());
  };
  const auto end = [](const SectionCurve& c) { return c.curve3d.evaluate(c.curve3d.domain_end()); };
  const auto reverse = [](SectionCurve& c) {
    c.curve3d = c.curve3d.reversed();
    c.pcurve = c.pcurve.reversed();
  };
  std::vector<bool> used(pieces.size(), false);
  for (std::size_t k = 0; k < pieces.size(); ++k) {
    if (used[k]) continue;
    used[k] = true;
    SectionLoop loop;
    loop.curves.push_back(pieces[k]);
    const Vec3 first = start(loop.curves.front());
    while (distance(end(loop.curves.back()), first) > opt.join_tolerance ||
           loop.curves.size() == 0) {
      const Vec3 e = end(loop.curves.back());
      std::size_t best = pieces.size();
      bool flip = false;
      double bd = opt.join_tolerance;
      for (std::size_t j = 0; j < pieces.size(); ++j) {
        if (used[j]) continue;
        const double ds = distance(start(pieces[j]), e), de = distance(end(pieces[j]), e);
        if (ds <= bd) {
          bd = ds;
          best = j;
          flip = false;
        }
        if (de <= bd) {
          bd = de;
          best = j;
          flip = true;
        }
      }
      if (best == pieces.size()) {
        char buf[160];
        std::snprintf(buf, sizeof buf,
                      "section: open chain at (%.9g, %.9g, %.9g) (open shell or tolerance)", e.x,
                      e.y, e.z);
        throw GeometryError(buf);
      }
      used[best] = true;
      SectionCurve c = pieces[best];
      if (flip) reverse(c);
      loop.curves.push_back(std::move(c));
    }
    // orientation: the solid lies on the left of the tangent (n x T points into the solid)
    {
      const SectionCurve& c = loop.curves.front();
      const double tm = 0.5 * (c.curve3d.domain_start() + c.curve3d.domain_end());
      const Vec3 T = c.curve3d.tangent(tm), N = c.normal(tm);
      if (dot(cross(n, T), N) > 0.0) {  // left side is outside
        std::reverse(loop.curves.begin(), loop.curves.end());
        for (auto& x : loop.curves) reverse(x);
      }
    }
    for (auto& x : loop.curves) x.curve = to_plane(x.curve3d, out.origin, out.s_axis, out.t_axis);
    loop.area = loop_area(loop.curves);
    out.loops.push_back(std::move(loop));
  }
  return out;
}

// ------------------------------------------------------------------ critical heights

namespace {

/// Bisection refinement of a sign change of g between a and b (g(a) has sign of ga).
template <class G>
double bisect(G&& g, double a, double b, double ga) {
  for (int it = 0; it < 200 && b - a > 1e-15 * std::max(1.0, std::abs(a)); ++it) {
    const double m = 0.5 * (a + b), gm = g(m);
    if ((gm < 0.0) == (ga < 0.0)) {
      a = m;
      ga = gm;
    } else {
      b = m;
    }
  }
  return 0.5 * (a + b);
}

/// Interior extrema of the height h = d . S(u, v) inside the trimmed face: for faces whose height
/// depends on one parameter the zeros of h'; otherwise Newton on (S_u . d, S_v . d) = 0 from seeds
/// on every knot span (Bezier patch).
void face_extrema(const TrimmedSurface& f, const Vec3& d, double size, std::vector<double>& out) {
  const NurbsSurface& s = f.surface();
  const auto du = s.domain_u(), dv = s.domain_v();
  for (const bool along_v : {true, false}) {
    const auto hc = height_curve(s, d, along_v, 1e-12 * size);
    if (!hc) continue;
    const auto kn = unique_knots(hc->knots()).first;
    const std::size_t n = 64 * std::max<std::size_t>(1, kn.size() - 1);
    const double t0 = hc->domain_start(), t1 = hc->domain_end();
    const auto dh = [&](double t) { return hc->derivatives(t, 1)[1].x; };
    double ta = t0, fa = dh(t0);
    for (std::size_t k = 1; k <= n; ++k) {
      const double tb = t0 + (t1 - t0) * static_cast<double>(k) / static_cast<double>(n),
                   fb = dh(tb);
      if ((fa < 0.0) != (fb < 0.0) && fa != 0.0 && fb != 0.0) {
        const double t = bisect(dh, ta, tb, fa);
        bool reached = f.is_untrimmed();
        for (int q = 0; q <= 64 && !reached; ++q) {
          const double w =
              along_v ? du[0] + (du[1] - du[0]) * q / 64.0 : dv[0] + (dv[1] - dv[0]) * q / 64.0;
          reached = along_v ? f.contains(w, t) : f.contains(t, w);
        }
        if (reached) out.push_back(hc->evaluate(t).x);
      }
      ta = tb;
      fa = fb;
    }
    return;
  }
  const auto ku = unique_knots(s.knots_u()).first, kv = unique_knots(s.knots_v()).first;
  for (std::size_t i = 0; i + 1 < ku.size(); ++i)
    for (std::size_t j = 0; j + 1 < kv.size(); ++j)
      for (int a = 0; a < 3; ++a)
        for (int b = 0; b < 3; ++b) {
          double u = ku[i] + (ku[i + 1] - ku[i]) * (a + 0.5) / 3.0,
                 v = kv[j] + (kv[j + 1] - kv[j]) * (b + 0.5) / 3.0;
          bool ok = false;
          for (int it = 0; it < 50; ++it) {
            const auto D = s.derivatives(u, v, 2);  // D[k][l]: d^(k+l) S / du^k dv^l
            const double gu = dot(D[1][0], d), gv = dot(D[0][1], d);
            const double huu = dot(D[2][0], d), huv = dot(D[1][1], d), hvv = dot(D[0][2], d);
            const double det = huu * hvv - huv * huv;
            if (std::abs(det) < 1e-300) break;
            const double su = (hvv * gu - huv * gv) / det, sv = (huu * gv - huv * gu) / det;
            u -= su;
            v -= sv;
            if (u < du[0] || u > du[1] || v < dv[0] || v > dv[1]) break;
            if (std::abs(su) + std::abs(sv) < 1e-14 * (1.0 + std::abs(u) + std::abs(v))) {
              ok = true;
              break;
            }
          }
          if (ok && f.contains(u, v)) out.push_back(dot(s.evaluate(u, v), d));
        }
}

}  // namespace

std::vector<double> critical_heights(const Shell& sh, const Vec3& dir, const Tolerance& tol) {
  const Vec3 d = normalized(dir);
  const double size = model_size(sh);
  std::vector<double> h;
  for (const auto& p : sh.vertices()) h.push_back(dot(p, d));
  for (const auto& f : sh.faces()) {
    const auto hv = height_curve(f.surface(), d, true, 1e-12 * size);
    const auto hu = height_curve(f.surface(), d, false, 1e-12 * size);
    if (hv && hu) {
      h.push_back(hv->control_points()[0].x);  // planar face normal to d
      continue;
    }
    face_extrema(f, d, size, h);
  }
  // extrema of the height along the edges: C' . d = 0
  for (const auto& e : sh.edges()) {
    const NurbsCurve& c = e.curve;
    const auto k = unique_knots(c.knots()).first;
    const std::size_t n = 64 * std::max<std::size_t>(1, k.size() - 1);
    const double t0 = c.domain_start(), t1 = c.domain_end();
    const auto g = [&](double t) { return dot(c.derivatives(t, 1)[1], d); };
    double ta = t0, fa = g(t0);
    for (std::size_t q = 1; q <= n; ++q) {
      const double tb = t0 + (t1 - t0) * static_cast<double>(q) / static_cast<double>(n),
                   fb = g(tb);
      if ((fa < 0.0) != (fb < 0.0) && fa != 0.0 && fb != 0.0)
        h.push_back(dot(c.evaluate(bisect(g, ta, tb, fa)), d));
      ta = tb;
      fa = fb;
    }
  }
  std::sort(h.begin(), h.end());
  std::vector<double> out;
  for (double x : h)
    if (out.empty() || x - out.back() > tol.length) out.push_back(x);
  return out;
}

bool is_extrusion_band(const Shell& sh, const Vec3& dir, double h0, double h1,
                       const Tolerance& tol) {
  const Vec3 d = normalized(dir);
  const double size = model_size(sh), margin = 1e-9 * size;
  for (const auto& f : sh.faces()) {
    const NurbsSurface& s = f.surface();
    // samples of the trimmed face strictly inside the band: normal perpendicular to d
    const auto du = s.domain_u(), dv = s.domain_v();
    const auto ku = unique_knots(s.knots_u()).first, kv = unique_knots(s.knots_v()).first;
    const std::size_t nu = 16 * std::max<std::size_t>(1, ku.size() - 1);
    const std::size_t nv = 16 * std::max<std::size_t>(1, kv.size() - 1);
    bool inBand = false;
    for (std::size_t i = 0; i <= nu; ++i)
      for (std::size_t j = 0; j <= nv; ++j) {
        const double u = du[0] + (du[1] - du[0]) * static_cast<double>(i) / static_cast<double>(nu);
        const double v = dv[0] + (dv[1] - dv[0]) * static_cast<double>(j) / static_cast<double>(nv);
        const double hh = dot(s.evaluate(u, v), d);
        if (hh <= h0 + margin || hh >= h1 - margin || !f.contains(u, v)) continue;
        inBand = true;
        if (std::abs(dot(s.normal(u, v), d)) > tol.angular) return false;
      }
    if (!inBand) continue;
    // control nets: every Bezier patch reaching into the band is a ruled extrusion along d (in one
    // parameter direction the control points differ by multiples of d only)
    for (const auto& b : s.decompose_bezier()) {
      double lo = 1e300, hi = -1e300;
      for (const auto& p : b.control_points()) {
        lo = std::min(lo, dot(p, d));
        hi = std::max(hi, dot(p, d));
      }
      if (hi <= h0 + margin || lo >= h1 - margin) continue;
      bool ruled = false;
      for (const bool along_u : {false, true}) {
        bool all = true;
        const std::size_t na = along_u ? b.num_v() : b.num_u(),
                          nb = along_u ? b.num_u() : b.num_v();
        for (std::size_t a = 0; a < na && all; ++a)
          for (std::size_t c = 1; c < nb && all; ++c) {
            const Vec3 p0 = along_u ? b.control_point(0, a) : b.control_point(a, 0);
            const Vec3 pc = along_u ? b.control_point(c, a) : b.control_point(a, c);
            const Vec3 diff = pc - p0;
            all = norm(diff - dot(diff, d) * d) <= 1e-9 * size;
          }
        ruled = ruled || all;
      }
      if (!ruled) return false;
    }
  }
  return true;
}

}  // namespace nf
