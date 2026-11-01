// Native STEP (ISO 10303-21) reader: Part 21 parser, geometry / topology translation and
// reconstruction of the trim loops in the parameter space of every face.
#include <algorithm>
#include <array>
#include <cctype>
#include <charconv>
#include <cmath>
#include <fstream>
#include <initializer_list>
#include <map>
#include <numbers>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

#include "nurbforge/analysis/closest_point.hpp"
#include "nurbforge/core/error.hpp"
#include "nurbforge/fitting/fit.hpp"
#include "nurbforge/geometry/primitives.hpp"
#include "nurbforge/geometry/transform.hpp"
#include "nurbforge/io/step.hpp"
#include "nurbforge/modeling/construct.hpp"
#include "nurbforge/nurbs/join.hpp"

namespace nf {

namespace {

constexpr double kPi = std::numbers::pi;
constexpr double kTwoPi = 2.0 * std::numbers::pi;

// ================================================================== Part 21 parser

struct Value {
  enum class Kind { Null, Derived, Number, String, Enum, Ref, List, Typed };
  Kind kind{Kind::Null};
  double num{0.0};
  std::size_t ref{0};
  std::string text;          // string contents, enumeration name or type keyword
  std::vector<Value> items;  // list items or typed parameters
};

struct Record {
  std::string name;
  std::vector<Value> params;
};

/// One record: simple entity instance; several: complex instance (one record per supertype).
using Instance = std::vector<Record>;
using InstanceMap = std::unordered_map<std::size_t, Instance>;

class Parser {
 public:
  explicit Parser(const std::string& text) : s_(text) {}

  InstanceMap parse() {
    InstanceMap out;
    bool data = false, started = false;
    for (;;) {
      skip();
      if (pos_ >= s_.size()) break;
      if (s_[pos_] == '#') {
        if (!data) fail("entity instance outside the DATA section");
        ++pos_;
        const std::size_t id = integer();
        skip();
        expect('=');
        skip();
        Instance inst;
        if (peek() == '(') {  // complex instance
          ++pos_;
          for (;;) {
            skip();
            if (peek() == ')') {
              ++pos_;
              break;
            }
            inst.push_back(record());
          }
          if (inst.empty()) fail("empty complex entity instance");
        } else {
          inst.push_back(record());
        }
        skip();
        expect(';');
        if (!out.emplace(id, std::move(inst)).second) fail("duplicate entity instance #" + std::to_string(id));
        continue;
      }
      const std::string kw = keyword();
      skip();
      if (kw == "ISO-10303-21") started = true;
      else if (!started) fail("not an ISO 10303-21 file (missing 'ISO-10303-21;')");
      if (kw == "DATA") {
        data = true;
        if (peek() == '(') (void)list();
      } else if (kw == "ENDSEC") {
        data = false;
      } else if (kw == "END-ISO-10303-21") {
        expect(';');
        break;
      } else if (kw != "ISO-10303-21" && kw != "HEADER") {
        (void)list();  // header entity (FILE_DESCRIPTION, FILE_NAME, FILE_SCHEMA, ...)
      }
      skip();
      expect(';');
    }
    return out;
  }

 private:
  const std::string& s_;
  std::size_t pos_{0};

  [[noreturn]] void fail(const std::string& msg) const {
    const auto line = 1 + std::count(s_.begin(), s_.begin() + static_cast<std::ptrdiff_t>(std::min(pos_, s_.size())), '\n');
    throw GeometryError("read_step: parse error at line " + std::to_string(line) + ": " + msg);
  }

  [[nodiscard]] char peek() const { return pos_ < s_.size() ? s_[pos_] : '\0'; }

  void expect(char c) {
    if (peek() != c) fail(std::string("expected '") + c + "'");
    ++pos_;
  }

  void skip() {
    while (pos_ < s_.size()) {
      const auto c = static_cast<unsigned char>(s_[pos_]);
      if (c <= ' ') {
        ++pos_;
      } else if (c == '/' && pos_ + 1 < s_.size() && s_[pos_ + 1] == '*') {
        const auto end = s_.find("*/", pos_ + 2);
        if (end == std::string::npos) fail("unterminated comment");
        pos_ = end + 2;
      } else {
        break;
      }
    }
  }

  static bool keyword_char(char c) {
    return std::isalnum(static_cast<unsigned char>(c)) != 0 || c == '_' || c == '-' || c == '!';
  }

  std::string keyword() {
    const std::size_t start = pos_;
    while (pos_ < s_.size() && keyword_char(s_[pos_])) ++pos_;
    if (pos_ == start) fail("expected a keyword");
    return s_.substr(start, pos_ - start);
  }

  std::size_t integer() {
    const std::size_t start = pos_;
    while (pos_ < s_.size() && std::isdigit(static_cast<unsigned char>(s_[pos_])) != 0) ++pos_;
    if (pos_ == start) fail("expected an instance number");
    std::size_t v = 0;
    std::from_chars(s_.data() + start, s_.data() + pos_, v);
    return v;
  }

  Record record() {
    Record r;
    r.name = keyword();
    skip();
    r.params = list();
    return r;
  }

  std::vector<Value> list() {
    expect('(');
    std::vector<Value> out;
    skip();
    if (peek() == ')') {
      ++pos_;
      return out;
    }
    for (;;) {
      skip();
      out.push_back(value());
      skip();
      if (peek() == ',') {
        ++pos_;
      } else {
        expect(')');
        return out;
      }
    }
  }

  Value value() {
    Value v;
    const char c = peek();
    if (c == '$') {
      ++pos_;
    } else if (c == '*') {
      v.kind = Value::Kind::Derived;
      ++pos_;
    } else if (c == '#') {
      ++pos_;
      v.kind = Value::Kind::Ref;
      v.ref = integer();
    } else if (c == '\'') {
      v.kind = Value::Kind::String;
      ++pos_;
      for (;;) {
        if (pos_ >= s_.size()) fail("unterminated string");
        if (s_[pos_] == '\'') {
          if (pos_ + 1 < s_.size() && s_[pos_ + 1] == '\'') {
            v.text += '\'';
            pos_ += 2;
            continue;
          }
          ++pos_;
          break;
        }
        v.text += s_[pos_++];
      }
    } else if (c == '"') {  // binary
      v.kind = Value::Kind::String;
      const auto end = s_.find('"', pos_ + 1);
      if (end == std::string::npos) fail("unterminated binary value");
      v.text = s_.substr(pos_ + 1, end - pos_ - 1);
      pos_ = end + 1;
    } else if (c == '.') {
      v.kind = Value::Kind::Enum;
      const auto end = s_.find('.', pos_ + 1);
      if (end == std::string::npos) fail("unterminated enumeration");
      v.text = s_.substr(pos_ + 1, end - pos_ - 1);
      pos_ = end + 1;
    } else if (c == '(') {
      v.kind = Value::Kind::List;
      v.items = list();
    } else if (std::isalpha(static_cast<unsigned char>(c)) != 0 || c == '!' || c == '_') {
      v.kind = Value::Kind::Typed;  // e.g. LENGTH_MEASURE(1.)
      v.text = keyword();
      skip();
      v.items = list();
    } else if (std::isdigit(static_cast<unsigned char>(c)) != 0 || c == '+' || c == '-') {
      v.kind = Value::Kind::Number;
      std::size_t start = pos_;
      while (pos_ < s_.size() && (std::isdigit(static_cast<unsigned char>(s_[pos_])) != 0 || s_[pos_] == '+' ||
                                  s_[pos_] == '-' || s_[pos_] == '.' || s_[pos_] == 'E' || s_[pos_] == 'e'))
        ++pos_;
      if (s_[start] == '+') ++start;
      const auto res = std::from_chars(s_.data() + start, s_.data() + pos_, v.num);
      if (res.ec != std::errc() || res.ptr != s_.data() + pos_) fail("invalid number");
    } else {
      fail(std::string("unexpected character '") + c + "'");
    }
    return v;
  }
};

// ================================================================== geometry helpers

using P2 = std::array<double, 2>;

struct Frame {
  Vec3 o, x, y, z;
};

double wrap_angle(double a) {  // [0, 2 pi)
  a = std::fmod(a, kTwoPi);
  if (a < 0.0) a += kTwoPi;
  return a >= kTwoPi ? 0.0 : a;
}

double wrap_signed(double a) {  // (-pi, pi]
  a = wrap_angle(a);
  return a > kPi ? a - kTwoPi : a;
}

P2 start2(const NurbsCurve& c) {
  const Vec3 p = c.evaluate(c.domain_start());
  return {p.x, p.y};
}
P2 end2(const NurbsCurve& c) {
  const Vec3 p = c.evaluate(c.domain_end());
  return {p.x, p.y};
}

NurbsCurve line2(const P2& a, const P2& b) { return make_line({a[0], a[1], 0.0}, {b[0], b[1], 0.0}); }

/// v / |v|, or the fallback for (nearly) zero vectors.
Vec3 unit(const Vec3& v, const Vec3& fallback) {
  const double n = norm(v);
  return n > 1e-300 ? v / n : fallback;
}

/// Exact rational quadratic arc of a hyperbola / parabola point(t) on [t0, t1] (tangent
/// direction deriv(t)): one Bezier segment per parameter interval of length <= 1; the middle
/// weight of a hyperbola segment is cosh of half its parameter span (1 for a parabola).
template <class Point, class Deriv>
NurbsCurve conic_arc(Point&& point, Deriv&& deriv, double t0, double t1, bool hyperbola) {
  const std::size_t n = hyperbola ? std::max<std::size_t>(1, static_cast<std::size_t>(std::ceil(t1 - t0))) : 1;
  std::vector<Vec3> P{point(t0)};
  std::vector<double> W{1.0};
  KnotVector K{0.0, 0.0, 0.0};
  for (std::size_t k = 0; k < n; ++k) {
    const double a = t0 + (t1 - t0) * static_cast<double>(k) / static_cast<double>(n);
    const double b = t0 + (t1 - t0) * static_cast<double>(k + 1) / static_cast<double>(n);
    P.push_back(detail::intersect_lines(point(a), deriv(a), point(b), deriv(b)));
    W.push_back(hyperbola ? std::cosh(0.5 * (b - a)) : 1.0);
    P.push_back(point(b));
    W.push_back(1.0);
    const double kv = static_cast<double>(k + 1);
    K.insert(K.end(), k + 1 < n ? std::size_t{2} : std::size_t{3}, kv);
  }
  return {2, std::move(K), std::move(P), std::move(W)};
}

/// Clamped curve with its first / last control point moved (end points set exactly).
NurbsCurve with_start(const NurbsCurve& c, const Vec3& p) {
  auto P = c.control_points();
  P.front() = p;
  return {c.degree(), c.knots(), std::move(P), c.weights()};
}
NurbsCurve with_end(const NurbsCurve& c, const Vec3& p) {
  auto P = c.control_points();
  P.back() = p;
  return {c.degree(), c.knots(), std::move(P), c.weights()};
}

std::size_t span_count(const NurbsCurve& c, double a, double b) {
  std::size_t n = 1;
  for (double k : unique_knots(c.knots()).first)
    if (k > a && k < b) ++n;
  return n;
}

/// Parameters on [a, b]: the curve's knots in (a, b) plus `per_span` intervals per span.
std::vector<double> sample_params(const NurbsCurve& c, double a, double b, std::size_t per_span) {
  std::vector<double> brk{a};
  for (double k : unique_knots(c.knots()).first)
    if (k > a + 1e-12 * (b - a) && k < b - 1e-12 * (b - a)) brk.push_back(k);
  brk.push_back(b);
  std::vector<double> t;
  for (std::size_t i = 0; i + 1 < brk.size(); ++i)
    for (std::size_t j = 0; j < per_span; ++j)
      t.push_back(brk[i] + (brk[i + 1] - brk[i]) * static_cast<double>(j) / static_cast<double>(per_span));
  t.push_back(b);
  return t;
}

/// Angular coordinate samples of a face boundary (for choosing the parameter patch of a surface
/// of revolution): angle, component of the inward boundary normal along the angle direction.
struct AngPt {
  double a{0.0}, s{0.0};
  bool ok{false};
};

/// True if a boundary polyline crosses the angle `mid`.
bool crosses(const std::vector<std::vector<AngPt>>& curves, double mid) {
  for (const auto& c : curves)
    for (std::size_t k = 0; k + 1 < c.size(); ++k) {
      if (!c[k].ok || !c[k + 1].ok) continue;
      const double d = wrap_signed(c[k + 1].a - c[k].a);
      if (std::abs(d) > 0.5 * kPi) continue;  // jump through a pole
      const double x = wrap_signed(mid - c[k].a);
      if ((d > 0.0 && x > 0.0 && x < d) || (d < 0.0 && x < 0.0 && x > d)) return true;
    }
  return false;
}

/// Angular range {start, sweep} covering the face: the complement of the largest angular gap
/// that lies outside the face (no boundary crosses it, the boundary turns away from it on both
/// sides). Full turn starting at `seam` (or `fallback`) if there is no such gap.
std::pair<double, double> angular_range(const std::vector<std::vector<AngPt>>& curves, std::optional<double> seam,
                                        double fallback) {
  std::vector<AngPt> all;
  for (const auto& c : curves)
    for (const auto& q : c)
      if (q.ok) all.push_back({wrap_angle(q.a), q.s, true});
  const double full_start = seam ? *seam : fallback;
  if (all.empty()) return {full_start, kTwoPi};
  std::sort(all.begin(), all.end(), [](const AngPt& p, const AngPt& q) { return p.a < q.a; });
  struct Group {
    double a, smin, smax;
  };
  std::vector<Group> g;
  for (const auto& q : all) {
    if (!g.empty() && q.a - g.back().a <= 1e-9) {
      g.back().smin = std::min(g.back().smin, q.s);
      g.back().smax = std::max(g.back().smax, q.s);
    } else {
      g.push_back({q.a, q.s, q.s});
    }
  }
  if (g.size() > 1 && g.front().a + kTwoPi - g.back().a <= 1e-9) {
    g.front().smin = std::min(g.front().smin, g.back().smin);
    g.front().smax = std::max(g.front().smax, g.back().smax);
    g.pop_back();
  }
  double best = -1.0, start = 0.0;
  for (std::size_t i = 0; i < g.size(); ++i) {
    const std::size_t j = (i + 1) % g.size();
    double gap = g[j].a - g[i].a;
    if (j <= i) gap += kTwoPi;
    if (gap < 1e-6 || gap <= best) continue;
    if (g[i].smax > 1e-6 || g[j].smin < -1e-6) continue;
    if (crosses(curves, g[i].a + 0.5 * gap)) continue;
    best = gap;
    start = g[j].a;
  }
  if (best < 0.0) return {full_start, kTwoPi};
  const double sweep = kTwoPi - best;
  const double margin = std::min(0.25 * best, 0.1);
  return {start - margin, sweep + 2.0 * margin};
}

// ================================================================== reader

struct Curve3 {
  enum class Kind { Line, Conic, Hyperbola, Parabola, Spline };
  Kind kind{Kind::Line};
  Vec3 p, d;               // line
  Frame f;                 // conic
  double a{0.0}, b{0.0};   // conic semi-axes (circle: a = b)
  std::optional<NurbsCurve> spline;
  bool reversed{false};    // TRIMMED_CURVE with sense_agreement .F.
};

struct Surf3 {
  enum class Kind { Plane, Cylinder, Cone, Sphere, Torus, Spline };
  Kind kind{Kind::Plane};
  Frame f;
  double r1{0.0}, r2{0.0};  // radius / major radius; cone: semi-angle, torus: minor radius
  std::optional<NurbsSurface> spline;
};

/// A face boundary loop in 3D: curves in loop direction (face on the left, seen against the
/// face normal) and the edge instances they come from.
struct Loop3 {
  std::vector<NurbsCurve> curves;
  std::vector<std::size_t> edges;
};

// ------------------------------------------------------------------ parameter-space mapping

/// Maps 3D boundary curves of a face into the parameter space of its surface.
class FaceMapper {
 public:
  FaceMapper(const NurbsSurface& s, bool face_left, double tol, std::vector<std::string>& notes)
      : S_(s), left_(face_left), tol_(tol), fit_tol_(0.1 * tol), notes_(notes) {
    const auto du = S_.domain_u(), dv = S_.domain_v();
    u0_ = du[0];
    u1_ = du[1];
    v0_ = dv[0];
    v1_ = dv[1];
    bool uc = true, vc = true;
    for (int k = 0; k <= 6; ++k) {
      const double fu = u0_ + (u1_ - u0_) * k / 6.0, fv = v0_ + (v1_ - v0_) * k / 6.0;
      uc = uc && distance(S_.evaluate(u0_, fv), S_.evaluate(u1_, fv)) <= 0.01 * tol_;
      vc = vc && distance(S_.evaluate(fu, v0_), S_.evaluate(fu, v1_)) <= 0.01 * tol_;
    }
    uc_ = uc;
    vc_ = vc;
    const auto pole_v = [&](double v) -> std::optional<Vec3> {
      const Vec3 p = S_.evaluate(0.5 * (u0_ + u1_), v);
      for (int k = 0; k <= 6; ++k)
        if (distance(S_.evaluate(u0_ + (u1_ - u0_) * k / 6.0, v), p) > 0.01 * tol_) return std::nullopt;
      return p;
    };
    const auto pole_u = [&](double u) -> std::optional<Vec3> {
      const Vec3 p = S_.evaluate(u, 0.5 * (v0_ + v1_));
      for (int k = 0; k <= 6; ++k)
        if (distance(S_.evaluate(u, v0_ + (v1_ - v0_) * k / 6.0), p) > 0.01 * tol_) return std::nullopt;
      return p;
    };
    pole_[0] = pole_u(u0_);
    pole_[1] = pole_u(u1_);
    pole_[2] = pole_v(v0_);
    pole_[3] = pole_v(v1_);
    // planar bilinear patch with affine parametrisation: exact inverse map
    if (S_.degree_u() == 1 && S_.degree_v() == 1 && S_.num_u() == 2 && S_.num_v() == 2) {
      const auto& W = S_.weights();
      const Vec3 p00 = S_.control_point(0, 0), p10 = S_.control_point(1, 0), p01 = S_.control_point(0, 1),
                 p11 = S_.control_point(1, 1);
      const Vec3 eu = p10 - p00, ev = p01 - p00;
      const double size = std::max(norm(eu), norm(ev));
      const bool equal_w = std::all_of(W.begin(), W.end(), [&](double w) { return std::abs(w - W[0]) <= 1e-14 * W[0]; });
      const double g11 = dot(eu, eu), g12 = dot(eu, ev), g22 = dot(ev, ev), det = g11 * g22 - g12 * g12;
      if (equal_w && norm(p00 + p11 - p10 - p01) <= 1e-13 * size && det > 1e-24 * g11 * g22) {
        const Vec3 ru = (u1_ - u0_) / det * (g22 * eu - g12 * ev);
        const Vec3 rv = (v1_ - v0_) / det * (g11 * ev - g12 * eu);
        affine_ = Transform({ru.x, ru.y, ru.z, u0_ - dot(ru, p00), rv.x, rv.y, rv.z, v0_ - dot(rv, p00), 0.0, 0.0,
                             0.0, 0.0});
      }
    }
  }

  /// The trimmed face bounded by the loops (untrimmed surface if there are none).
  TrimmedSurface build(const std::vector<Loop3>& loops) {
    if (loops.empty()) return TrimmedSurface(S_);
    std::vector<TrimLoop> closed;
    std::vector<std::pair<TrimLoop, int>> open;  // loops spanning a period (0: in u, 1: in v)
    for (const auto& l : loops) {
      std::vector<NurbsCurve> pcs;
      for (const auto& c : l.curves)
        for (auto& pc : map_curve(c)) pcs.push_back(std::move(pc));
      if (pcs.empty()) continue;
      const std::size_t m = pcs.size();
      std::vector<std::optional<NurbsCurve>> insert(m);
      std::size_t period_at = m;
      int period_dir = 0, periods = 0;
      for (std::size_t i = 0; i < m; ++i) {
        const std::size_t j = (i + 1) % m;
        const P2 e = end2(pcs[i]), s = start2(pcs[j]);
        switch (const int kind = junction(e, s); kind) {
          case 0: {  // snap
            const Vec3 mid{0.5 * (e[0] + s[0]), 0.5 * (e[1] + s[1]), 0.0};
            pcs[i] = with_end(pcs[i], mid);
            pcs[j] = with_start(pcs[j], mid);
            break;
          }
          case 1:  // degenerate in 3D (pole): straight connection
            insert[i] = line2(e, s);
            break;
          default:  // across the seam of a closed surface
            ++periods;
            period_at = i;
            period_dir = kind - 2;
            break;
        }
      }
      if (periods > 1) throw GeometryError("face boundary crosses the seam of a closed surface several times");
      TrimLoop tl;
      const std::size_t first = periods ? (period_at + 1) % m : 0;
      for (std::size_t k = 0; k < m; ++k) {
        const std::size_t i = (first + k) % m;
        tl.curves.push_back(pcs[i]);
        if (insert[i]) tl.curves.push_back(*insert[i]);
      }
      if (periods) open.emplace_back(std::move(tl), period_dir);
      else closed.push_back(std::move(tl));
    }
    close_open_loops(open, closed);
    if (closed.empty()) return TrimmedSurface(S_);
    std::size_t outer = 0;
    double best = -1.0;
    for (std::size_t i = 0; i < closed.size(); ++i) {
      const double a = std::abs(closed[i].signed_area());
      if (a > best) {
        best = a;
        outer = i;
      }
    }
    std::vector<TrimLoop> holes;
    for (std::size_t i = 0; i < closed.size(); ++i)
      if (i != outer) holes.push_back(std::move(closed[i]));
    return {S_, std::move(closed[outer]), std::move(holes)};
  }

 private:
  struct Smp {
    double t{0.0};
    double u{0.0}, v{0.0}, d{0.0};
    bool amb_u{false}, amb_v{false}, free_u{false}, free_v{false};
  };

  const NurbsSurface& S_;
  bool left_;
  double tol_, fit_tol_;
  std::vector<std::string>& notes_;
  double u0_{0.0}, u1_{1.0}, v0_{0.0}, v1_{1.0};
  bool uc_{false}, vc_{false};
  std::array<std::optional<Vec3>, 4> pole_;  // at u0, u1, v0, v1
  std::optional<Transform> affine_;

  [[nodiscard]] Vec3 at(const P2& p) const { return S_.evaluate(p[0], p[1]); }

  /// Junction classification: 0 snap, 1 degenerate connection, 2 / 3 period jump in u / v.
  int junction(const P2& e, const P2& s) const {
    const double du = s[0] - e[0], dv = s[1] - e[1], gap = std::hypot(du, dv);
    const double scale = std::max(u1_ - u0_, v1_ - v0_);
    if (gap <= 1e-13 * scale) return 0;
    const Vec3 pe = at(e);
    bool degenerate = true;
    for (int k = 1; k <= 8 && degenerate; ++k)
      degenerate = distance(at({e[0] + du * k / 8.0, e[1] + dv * k / 8.0}), pe) <= tol_;
    if (degenerate) return gap <= 1e-9 * scale ? 0 : 1;
    const bool same3d = distance(pe, at(s)) <= tol_;
    if (uc_ && same3d && std::abs(std::abs(du) - (u1_ - u0_)) <= 1e-6 * (u1_ - u0_) && std::abs(dv) <= 1e-6 * (v1_ - v0_))
      return 2;
    if (vc_ && same3d && std::abs(std::abs(dv) - (v1_ - v0_)) <= 1e-6 * (v1_ - v0_) && std::abs(du) <= 1e-6 * (u1_ - u0_))
      return 3;
    throw GeometryError("trim loop is not closed in the parameter space of the surface (gap " +
                        std::to_string(distance(pe, at(s))) + ")");
  }

  /// Loops running once around a closed direction: one is closed through a pole, two are joined
  /// along the seam into one loop.
  void close_open_loops(std::vector<std::pair<TrimLoop, int>>& open, std::vector<TrimLoop>& closed) const {
    if (open.empty()) return;
    if (open.size() > 2 || (open.size() == 2 && open[0].second != open[1].second))
      throw GeometryError("unsupported face topology on a closed surface (" + std::to_string(open.size()) +
                          " loops around the seam)");
    const auto ends = [](const TrimLoop& l) { return std::make_pair(start2(l.curves.front()), end2(l.curves.back())); };
    if (open.size() == 1) {
      const TrimLoop& l = open[0].first;
      const auto [s, e] = ends(l);
      TrimLoop out = l;
      if (open[0].second == 0) {  // around u: close through the pole on the face side
        const bool plus_v = (e[0] > s[0]) == left_;
        const std::optional<Vec3>& pole = pole_[plus_v ? 3 : 2];
        if (!pole) throw GeometryError("face boundary on a closed surface does not enclose a region");
        const double vp = plus_v ? v1_ : v0_;
        out.curves.push_back(line2(e, {e[0], vp}));
        out.curves.push_back(line2({e[0], vp}, {s[0], vp}));
        out.curves.push_back(line2({s[0], vp}, s));
      } else {
        const bool minus_u = (e[1] > s[1]) == left_;
        const std::optional<Vec3>& pole = pole_[minus_u ? 0 : 1];
        if (!pole) throw GeometryError("face boundary on a closed surface does not enclose a region");
        const double up = minus_u ? u0_ : u1_;
        out.curves.push_back(line2(e, {up, e[1]}));
        out.curves.push_back(line2({up, e[1]}, {up, s[1]}));
        out.curves.push_back(line2({up, s[1]}, s));
      }
      closed.push_back(std::move(out));
      return;
    }
    TrimLoop a = open[0].first, b = open[1].first;
    const auto dir = static_cast<std::size_t>(open[0].second);
    const auto [as, ae] = ends(a);
    const auto [bs0, be0] = ends(b);
    if ((ae[dir] - as[dir] > 0.0) == (be0[dir] - bs0[dir] > 0.0)) b = b.reversed();
    const auto [bs, be] = ends(b);
    TrimLoop out = a;
    out.curves.push_back(line2(ae, bs));
    for (const auto& c : b.curves) out.curves.push_back(c);
    out.curves.push_back(line2(be, as));
    closed.push_back(std::move(out));
  }

  // ---------------------------------------------------------------- curve mapping

  std::vector<NurbsCurve> map_curve(const NurbsCurve& c) {
    if (affine_) {
      auto pc = c.transformed(*affine_);
      return {pc};
    }
    return map_range(c, c.domain_start(), c.domain_end(), 0);
  }

  SurfacePoint project(const Vec3& p, const Smp* prev) const {
    if (prev) {
      const SurfacePoint loc = closest_point_local(S_, p, prev->u, prev->v);
      if (loc.distance <= 10.0 * tol_) return loc;
      const SurfacePoint glob = closest_point(S_, p);
      return glob.distance <= loc.distance ? glob : loc;
    }
    return closest_point(S_, p);
  }

  std::vector<Smp> sample(const NurbsCurve& c, const std::vector<double>& ts) const {
    std::vector<Smp> out;
    out.reserve(ts.size());
    for (double t : ts) {
      const Vec3 p = c.evaluate(t);
      Smp s;
      s.t = t;
      s.u = 0.5 * (u0_ + u1_);
      s.v = 0.5 * (v0_ + v1_);
      // points at a pole are not projected (the free coordinate is resolved below)
      for (int k = 0; k < 4; ++k) {
        const auto& pole = pole_[static_cast<std::size_t>(k)];
        if (!pole || distance(p, *pole) > tol_) continue;
        s.d = distance(p, *pole);
        if (k < 2) {
          s.u = k == 0 ? u0_ : u1_;
          s.free_v = true;
        } else {
          s.v = k == 2 ? v0_ : v1_;
          s.free_u = true;
        }
      }
      if (!s.free_u && !s.free_v) {
        const Smp* prev = nullptr;
        for (auto it = out.rbegin(); it != out.rend() && !prev; ++it)
          if (!it->free_u && !it->free_v) prev = &*it;
        const SurfacePoint sp = project(p, prev);
        s.u = sp.u;
        s.v = sp.v;
        s.d = sp.distance;
      }
      const Vec3 q = S_.evaluate(s.u, s.v);
      const double eps = 0.01 * fit_tol_;
      s.amb_u = uc_ && !s.free_u &&
                (distance(S_.evaluate(u0_, s.v), q) <= eps || distance(S_.evaluate(u1_, s.v), q) <= eps);
      s.amb_v = vc_ && !s.free_v &&
                (distance(S_.evaluate(s.u, v0_), q) <= eps || distance(S_.evaluate(s.u, v1_), q) <= eps);
      out.push_back(s);
    }
    resolve(out);
    return out;
  }

  /// Seam and pole ambiguities: seam samples take the side of their nearest unambiguous
  /// neighbour (curves lying on the seam: the side given by the face orientation: the face lies
  /// to the left of the boundary iff left_), pole samples the free coordinate of their nearest
  /// regular neighbour.
  void resolve(std::vector<Smp>& smp) const {
    const std::size_t n = smp.size();
    const auto nearest = [&](std::size_t i, auto&& good) -> std::optional<std::size_t> {
      for (std::size_t d = 1; d < n; ++d) {
        if (i >= d && good(smp[i - d])) return i - d;
        if (i + d < n && good(smp[i + d])) return i + d;
      }
      return std::nullopt;
    };
    const auto pick = [](double x, double lo, double hi) { return std::abs(x - lo) < std::abs(x - hi) ? lo : hi; };
    std::vector<std::size_t> rest_u, rest_v;
    for (std::size_t i = 0; i < n; ++i) {
      if (smp[i].amb_u) {
        if (const auto j = nearest(i, [](const Smp& s) { return !s.amb_u && !s.free_u; }))
          smp[i].u = pick(smp[*j].u, u0_, u1_);
        else rest_u.push_back(i);
      }
      if (smp[i].amb_v) {
        if (const auto j = nearest(i, [](const Smp& s) { return !s.amb_v && !s.free_v; }))
          smp[i].v = pick(smp[*j].v, v0_, v1_);
        else rest_v.push_back(i);
      }
    }
    if (!rest_u.empty()) {  // moving +v with the face on the left: face at smaller u, seam u1
      const double dv = smp.back().v - smp.front().v;
      const double seam = (dv > 0.0) == left_ ? u1_ : u0_;
      for (std::size_t i : rest_u) smp[i].u = seam;
    }
    if (!rest_v.empty()) {  // moving +u with the face on the left: face at larger v, seam v0
      const double du = smp.back().u - smp.front().u;
      const double seam = (du > 0.0) == left_ ? v0_ : v1_;
      for (std::size_t i : rest_v) smp[i].v = seam;
    }
    for (std::size_t i = 0; i < n; ++i) {
      if (smp[i].free_u) {
        if (const auto j = nearest(i, [](const Smp& s) { return !s.free_u; })) smp[i].u = smp[*j].u;
        else throw GeometryError("degenerate edge at a pole");
      }
      if (smp[i].free_v) {
        if (const auto j = nearest(i, [](const Smp& s) { return !s.free_v; })) smp[i].v = smp[*j].v;
        else throw GeometryError("degenerate edge at a pole");
      }
    }
  }

  /// Parameter at which the curve must be split: passage through a pole or across a seam.
  std::optional<double> split_point(const NurbsCurve& c, const std::vector<Smp>& smp) const {
    const double a = smp.front().t, b = smp.back().t, margin = 1e-9 * (b - a);
    for (std::size_t k = 1; k + 1 < smp.size(); ++k)
      if (smp[k].free_u || smp[k].free_v) return smp[k].t;
    for (std::size_t k = 0; k + 1 < smp.size(); ++k) {
      const double du = std::abs(smp[k + 1].u - smp[k].u), dv = std::abs(smp[k + 1].v - smp[k].v);
      if (du <= 0.25 * (u1_ - u0_) && dv <= 0.25 * (v1_ - v0_)) continue;
      const double ta = smp[k].t, tb = smp[k + 1].t;
      for (const auto& pole : pole_) {  // closest approach to the pole (golden section)
        if (!pole) continue;
        double lo = ta, hi = tb;
        const double g = 0.5 * (std::sqrt(5.0) - 1.0);
        for (int it = 0; it < 80; ++it) {
          const double m1 = hi - g * (hi - lo), m2 = lo + g * (hi - lo);
          if (distance(c.evaluate(m1), *pole) < distance(c.evaluate(m2), *pole)) hi = m2;
          else lo = m1;
        }
        const double t = 0.5 * (lo + hi);
        if (distance(c.evaluate(t), *pole) <= tol_ && t > a + margin && t < b - margin) return t;
      }
      const bool seam_u = uc_ && du > 0.5 * (u1_ - u0_), seam_v = vc_ && dv > 0.5 * (v1_ - v0_);
      if (!seam_u && !seam_v) continue;
      const auto side = [&](double t) {
        const SurfacePoint sp = closest_point(S_, c.evaluate(t));
        return seam_u ? sp.u > 0.5 * (u0_ + u1_) : sp.v > 0.5 * (v0_ + v1_);
      };
      const bool side_a = seam_u ? smp[k].u > 0.5 * (u0_ + u1_) : smp[k].v > 0.5 * (v0_ + v1_);
      double lo = ta, hi = tb;
      for (int it = 0; it < 60 && hi - lo > 1e-13 * (b - a); ++it) {
        const double m = 0.5 * (lo + hi);
        if (side(m) == side_a) lo = m;
        else hi = m;
      }
      const double t = 0.5 * (lo + hi);
      if (t > a + margin && t < b - margin) return t;
    }
    return std::nullopt;
  }

  struct Fit {
    NurbsCurve curve;
    bool ok;
  };

  Fit fit(const NurbsCurve& c, const std::vector<Smp>& smp) const {
    std::vector<Vec3> pts;
    std::vector<double> ts;
    for (const auto& s : smp) {
      pts.push_back({s.u, s.v, 0.0});
      ts.push_back(s.t);
    }
    const std::size_t degree = std::min<std::size_t>(3, pts.size() - 1);
    NurbsCurve pc = interpolate_curve(pts, degree, ts);
    auto P = pc.control_points();
    for (auto& p : P) p = {std::clamp(p.x, u0_, u1_), std::clamp(p.y, v0_, v1_), 0.0};
    pc = NurbsCurve(pc.degree(), pc.knots(), std::move(P), pc.weights());
    bool ok = true;
    for (std::size_t k = 0; k + 1 < smp.size() && ok; ++k) {
      const double tm = 0.5 * (smp[k].t + smp[k + 1].t);
      const Vec3 q = pc.evaluate(tm);
      const double err = closest_point_local(c, S_.evaluate(q.x, q.y), tm).distance;
      ok = err <= fit_tol_ + 1.01 * std::max(smp[k].d, smp[k + 1].d);
    }
    return {std::move(pc), ok};
  }

  std::vector<NurbsCurve> map_range(const NurbsCurve& c, double a, double b, int depth) {
    const std::size_t spans = span_count(c, a, b);
    std::size_t per_span = std::max<std::size_t>(2, std::min<std::size_t>(16, 32 / spans));
    for (int attempt = 0;; ++attempt) {
      const auto smp = sample(c, sample_params(c, a, b, per_span));
      if (depth < 8)
        if (const auto t = split_point(c, smp)) {
          auto out = map_range(c, a, *t, depth + 1);
          for (auto& pc : map_range(c, *t, b, depth + 1)) out.push_back(std::move(pc));
          return out;
        }
      Fit f = fit(c, smp);
      if (f.ok) return {std::move(f.curve)};
      if (attempt >= 5 || smp.size() > 2048) {
        notes_.push_back("trim curve not within tolerance after refinement");
        return {std::move(f.curve)};
      }
      per_span *= 2;
    }
  }
};

// ------------------------------------------------------------------ entity translation

class Reader {
 public:
  Reader(InstanceMap inst, const StepReadOptions& opt) : inst_(std::move(inst)), opt_(opt), tol_(opt.tolerance) {
    if (!(tol_ > 0.0)) throw GeometryError("read_step: tolerance must be positive");
    read_units();
  }

  std::vector<Shell> read() {
    std::vector<std::size_t> ids;
    for (const auto& [id, in] : inst_) ids.push_back(id);
    std::sort(ids.begin(), ids.end());
    std::vector<Shell> out;
    for (std::size_t id : ids) {
      const Instance& in = inst_.at(id);
      if (in.size() != 1) continue;
      const Record& r = in[0];
      if (r.name == "MANIFOLD_SOLID_BREP" || r.name == "BREP_WITH_VOIDS") {
        std::vector<TrimmedSurface> faces;
        shell_faces(ref(arg(id, r, 1), id), false, faces);
        if (r.name == "BREP_WITH_VOIDS")
          for (const auto& v : items(arg(id, r, 2), id)) shell_faces(ref(v, id), false, faces);
        add_shell(out, std::move(faces));
      } else if (r.name == "SHELL_BASED_SURFACE_MODEL") {
        for (const auto& s : items(arg(id, r, 1), id)) {
          std::vector<TrimmedSurface> faces;
          shell_faces(ref(s, id), false, faces);
          add_shell(out, std::move(faces));
        }
      } else if (r.name == "FACETED_BREP" || r.name == "FACETED_BREP_SHAPE_REPRESENTATION") {
        problem("unsupported entity " + r.name + " (#" + std::to_string(id) + ")");
      }
    }
    return out;
  }

 private:
  InstanceMap inst_;
  const StepReadOptions& opt_;
  double tol_;
  double scale_{1.0};        // file length unit -> target unit
  double angle_unit_{1.0};   // file plane angle unit -> radians
  std::unordered_map<std::size_t, Vec3> points_;
  std::unordered_map<std::size_t, NurbsCurve> edges_;

  // ---------------------------------------------------------------- access helpers

  void problem(const std::string& msg) {
    if (!opt_.skip_unsupported) throw GeometryError("read_step: " + msg);
    if (opt_.warnings) opt_.warnings->push_back("skipped: " + msg);
  }

  const Instance& inst(std::size_t id) const {
    const auto it = inst_.find(id);
    if (it == inst_.end()) throw GeometryError("reference to missing entity instance #" + std::to_string(id));
    return it->second;
  }

  std::string type_of(std::size_t id) const {
    const Instance& in = inst(id);
    if (in.size() == 1) return in[0].name;
    std::string s = "(";
    for (std::size_t i = 0; i < in.size(); ++i) s += (i ? " " : "") + in[i].name;
    return s + ")";
  }

  [[noreturn]] void unsupported(std::size_t id) const {
    throw GeometryError("unsupported entity " + type_of(id) + " (#" + std::to_string(id) + ")");
  }

  [[noreturn]] void bad(std::size_t id, const std::string& msg) const {
    throw GeometryError(msg + " in " + type_of(id) + " (#" + std::to_string(id) + ")");
  }

  static const Record* find(const Instance& in, std::string_view name) {
    for (const auto& r : in)
      if (r.name == name) return &r;
    return nullptr;
  }

  /// The record of a simple instance (throws for other types).
  const Record& simple(std::size_t id, std::initializer_list<std::string_view> names) const {
    const Instance& in = inst(id);
    if (in.size() == 1)
      for (auto n : names)
        if (in[0].name == n) return in[0];
    unsupported(id);
  }

  const Value& arg(std::size_t id, const Record& r, std::size_t i) const {
    if (i >= r.params.size()) bad(id, "missing attribute " + std::to_string(i + 1));
    return r.params[i];
  }

  double num(const Value& v, std::size_t id) const {
    if (v.kind == Value::Kind::Number) return v.num;
    if (v.kind == Value::Kind::Typed && v.items.size() == 1) return num(v.items[0], id);
    bad(id, "expected a number");
  }

  std::size_t ref(const Value& v, std::size_t id) const {
    if (v.kind != Value::Kind::Ref) bad(id, "expected an entity reference");
    return v.ref;
  }

  bool flag(const Value& v, std::size_t id) const {
    if (v.kind == Value::Kind::Enum) {
      if (v.text == "T" || v.text == "U") return true;
      if (v.text == "F") return false;
    }
    bad(id, "expected a logical value");
  }

  const std::vector<Value>& items(const Value& v, std::size_t id) const {
    if (v.kind != Value::Kind::List) bad(id, "expected a list");
    return v.items;
  }

  std::size_t count(const Value& v, std::size_t id) const {
    const double x = num(v, id);
    if (!(x >= 0.0) || x != std::floor(x) || x > 1e9) bad(id, "expected a non-negative integer");
    return static_cast<std::size_t>(x);
  }

  // ---------------------------------------------------------------- units

  /// Size of a unit in SI base units (metre / radian).
  double unit_factor(std::size_t id, int depth = 0) const {
    if (depth > 8) bad(id, "recursive unit definition");
    const Instance& in = inst(id);
    if (const Record* si = find(in, "SI_UNIT")) {
      const std::size_t o = si->params.size() >= 3 ? 1 : 0;  // simple instance: dimensions first
      if (si->params.size() < o + 2) bad(id, "malformed SI_UNIT");
      double f = 1.0;
      const Value& prefix = si->params[o];
      if (prefix.kind == Value::Kind::Enum) {
        static const std::map<std::string, double> prefixes{
            {"EXA", 1e18},  {"PETA", 1e15}, {"TERA", 1e12},  {"GIGA", 1e9},   {"MEGA", 1e6},   {"KILO", 1e3},
            {"HECTO", 1e2}, {"DECA", 1e1},  {"DECI", 1e-1},  {"CENTI", 1e-2}, {"MILLI", 1e-3}, {"MICRO", 1e-6},
            {"NANO", 1e-9}, {"PICO", 1e-12}, {"FEMTO", 1e-15}, {"ATTO", 1e-18}};
        const auto it = prefixes.find(prefix.text);
        if (it == prefixes.end()) bad(id, "unknown SI prefix ." + prefix.text + ".");
        f = it->second;
      }
      return f;
    }
    if (const Record* cb = find(in, "CONVERSION_BASED_UNIT")) {
      if (cb->params.size() < 2) bad(id, "malformed CONVERSION_BASED_UNIT");
      const std::size_t m = ref(cb->params[1], id);
      const Instance& mi = inst(m);
      const Record* mr = nullptr;
      for (const auto& r : mi)
        if (r.params.size() >= 2 && r.params[1].kind == Value::Kind::Ref) mr = &r;
      if (!mr) bad(m, "malformed measure with unit");
      return num(mr->params[0], m) * unit_factor(mr->params[1].ref, depth + 1);
    }
    unsupported(id);
  }

  void read_units() {
    std::optional<double> length, angle;
    std::vector<std::size_t> ids;
    for (const auto& [id, in] : inst_)
      if (find(in, "GLOBAL_UNIT_ASSIGNED_CONTEXT")) ids.push_back(id);
    std::sort(ids.begin(), ids.end());
    for (std::size_t id : ids) {
      const Record* r = find(inst_.at(id), "GLOBAL_UNIT_ASSIGNED_CONTEXT");
      if (r->params.empty()) continue;
      for (const auto& u : items(r->params[0], id)) {
        const std::size_t uid = ref(u, id);
        const Instance& ui = inst(uid);
        if (!length && find(ui, "LENGTH_UNIT")) length = unit_factor(uid);
        if (!angle && find(ui, "PLANE_ANGLE_UNIT")) angle = unit_factor(uid);
      }
    }
    const double target = [&] {
      switch (opt_.unit) {
        case StepLengthUnit::Millimetre: return 1e-3;
        case StepLengthUnit::Centimetre: return 1e-2;
        case StepLengthUnit::Metre: return 1.0;
        case StepLengthUnit::Inch: return 0.0254;
      }
      return 1e-3;
    }();
    scale_ = length.value_or(1e-3) / target;
    angle_unit_ = angle.value_or(1.0);
  }

  double len(const Value& v, std::size_t id) const { return num(v, id) * scale_; }

  // ---------------------------------------------------------------- basic geometry

  Vec3 point(std::size_t id) {
    if (const auto it = points_.find(id); it != points_.end()) return it->second;
    const Record& r = simple(id, {"CARTESIAN_POINT"});
    const auto& c = items(arg(id, r, 1), id);
    if (c.size() != 3) bad(id, "expected a 3D point");
    const Vec3 p{num(c[0], id) * scale_, num(c[1], id) * scale_, num(c[2], id) * scale_};
    points_.emplace(id, p);
    return p;
  }

  Vec3 direction(std::size_t id) const {
    const Record& r = simple(id, {"DIRECTION"});
    const auto& c = items(arg(id, r, 1), id);
    if (c.size() != 3) bad(id, "expected a 3D direction");
    const Vec3 d{num(c[0], id), num(c[1], id), num(c[2], id)};
    if (!(norm(d) > 0.0)) bad(id, "zero direction");
    return d / norm(d);
  }

  Frame placement(std::size_t id) {
    const Record& r = simple(id, {"AXIS2_PLACEMENT_3D"});
    Frame f;
    f.o = point(ref(arg(id, r, 1), id));
    f.z = r.params.size() > 2 && r.params[2].kind == Value::Kind::Ref ? direction(r.params[2].ref) : Vec3{0, 0, 1};
    Vec3 x = r.params.size() > 3 && r.params[3].kind == Value::Kind::Ref ? direction(r.params[3].ref)
                                                                        : any_perpendicular(f.z);
    x = x - dot(x, f.z) * f.z;
    if (norm(x) < 1e-9) x = any_perpendicular(f.z);
    f.x = unit(x, any_perpendicular(f.z));
    f.y = cross(f.z, f.x);
    return f;
  }

  Vec3 vertex(std::size_t id) {
    const Record& r = simple(id, {"VERTEX_POINT"});
    return point(ref(arg(id, r, 1), id));
  }

  std::vector<Vec3> points(const Value& v, std::size_t id) {
    std::vector<Vec3> out;
    for (const auto& p : items(v, id)) out.push_back(point(ref(p, id)));
    return out;
  }

  /// Knot vector of a B-spline (explicit knots or implied by the curve / surface subtype).
  KnotVector knots(std::size_t id, const std::string& kind, std::size_t p, std::size_t n, const Value* mults,
                   const Value* values) const {
    KnotVector k;
    if (kind == "WITH_KNOTS") {
      const auto& m = items(*mults, id);
      const auto& v = items(*values, id);
      if (m.size() != v.size()) bad(id, "knot multiplicities and values differ in length");
      for (std::size_t i = 0; i < m.size(); ++i) k.insert(k.end(), count(m[i], id), num(v[i], id));
    } else if (kind == "BEZIER") {
      if (p == 0 || (n - 1) % p != 0) bad(id, "invalid number of Bezier control points");
      const std::size_t segs = (n - 1) / p;
      k.insert(k.end(), p + 1, 0.0);
      for (std::size_t s = 1; s < segs; ++s) k.insert(k.end(), p, static_cast<double>(s));
      k.insert(k.end(), p + 1, static_cast<double>(segs));
    } else if (kind == "QUASI_UNIFORM") {
      k = clamped_uniform_knots(p, n);
    } else if (kind == "UNIFORM") {
      for (std::size_t i = 0; i < n + p + 1; ++i) k.push_back(static_cast<double>(i) - static_cast<double>(p));
    } else {
      unsupported(id);
    }
    if (k.size() != n + p + 1) bad(id, "inconsistent knot vector");
    return k;
  }

  /// B-spline subtype of a curve / surface instance (prefix "B_SPLINE_CURVE" or "B_SPLINE_SURFACE").
  static std::string bspline_kind(const Instance& in, const std::string& base) {
    const std::string stem = base.substr(std::string("B_SPLINE_").size());  // CURVE / SURFACE
    for (const auto& r : in) {
      if (r.name == base + "_WITH_KNOTS") return "WITH_KNOTS";
      if (r.name == "BEZIER_" + stem) return "BEZIER";
      if (r.name == "QUASI_UNIFORM_" + stem) return "QUASI_UNIFORM";
      if (r.name == "UNIFORM_" + stem) return "UNIFORM";
    }
    return {};
  }

  NurbsCurve bspline_curve(std::size_t id) {
    const Instance& in = inst(id);
    const std::string kind = bspline_kind(in, "B_SPLINE_CURVE");
    if (kind.empty()) unsupported(id);
    std::vector<const Value*> base, kn;
    const Value* w = nullptr;
    if (in.size() == 1) {  // simple: name, 5 B_SPLINE_CURVE attributes, subtype attributes
      const Record& r = in[0];
      if (r.params.size() < 6) bad(id, "too few attributes");
      for (std::size_t i = 1; i < r.params.size(); ++i) (i < 6 ? base : kn).push_back(&r.params[i]);
    } else {
      const Record* b = find(in, "B_SPLINE_CURVE");
      if (!b || b->params.size() < 5) bad(id, "missing B_SPLINE_CURVE attributes");
      for (const auto& v : b->params) base.push_back(&v);
      if (const Record* k = find(in, "B_SPLINE_CURVE_WITH_KNOTS"))
        for (const auto& v : k->params) kn.push_back(&v);
      if (const Record* r = find(in, "RATIONAL_B_SPLINE_CURVE")) {
        if (r->params.empty()) bad(id, "missing weights");
        w = &r->params[0];
      }
    }
    if (kind == "WITH_KNOTS" && kn.size() < 2) bad(id, "missing knots");
    const std::size_t p = count(*base[0], id);
    std::vector<Vec3> P = points(*base[1], id);
    KnotVector K = knots(id, kind, p, P.size(), kind == "WITH_KNOTS" ? kn[0] : nullptr,
                         kind == "WITH_KNOTS" ? kn[1] : nullptr);
    std::vector<double> W;
    if (w)
      for (const auto& x : items(*w, id)) W.push_back(num(x, id));
    if (!W.empty() && W.size() != P.size()) bad(id, "number of weights differs from control points");
    return NurbsCurve(p, std::move(K), std::move(P), std::move(W)).clamped();
  }

  NurbsSurface bspline_surface(std::size_t id) {
    const Instance& in = inst(id);
    const std::string kind = bspline_kind(in, "B_SPLINE_SURFACE");
    if (kind.empty()) unsupported(id);
    std::vector<const Value*> base, kn;
    const Value* w = nullptr;
    if (in.size() == 1) {  // simple: name, 7 B_SPLINE_SURFACE attributes, subtype attributes
      const Record& r = in[0];
      if (r.params.size() < 8) bad(id, "too few attributes");
      for (std::size_t i = 1; i < r.params.size(); ++i) (i < 8 ? base : kn).push_back(&r.params[i]);
    } else {
      const Record* b = find(in, "B_SPLINE_SURFACE");
      if (!b || b->params.size() < 7) bad(id, "missing B_SPLINE_SURFACE attributes");
      for (const auto& v : b->params) base.push_back(&v);
      if (const Record* k = find(in, "B_SPLINE_SURFACE_WITH_KNOTS"))
        for (const auto& v : k->params) kn.push_back(&v);
      if (const Record* r = find(in, "RATIONAL_B_SPLINE_SURFACE")) {
        if (r->params.empty()) bad(id, "missing weights");
        w = &r->params[0];
      }
    }
    if (kind == "WITH_KNOTS" && kn.size() < 4) bad(id, "missing knots");
    const std::size_t pu = count(*base[0], id), pv = count(*base[1], id);
    const auto& rows = items(*base[2], id);
    if (rows.empty()) bad(id, "empty control net");
    const std::size_t nu = rows.size(), nv = items(rows[0], id).size();
    std::vector<Vec3> P;
    for (const auto& row : rows) {
      auto pts = points(row, id);
      if (pts.size() != nv) bad(id, "ragged control net");
      P.insert(P.end(), pts.begin(), pts.end());
    }
    const bool wk = kind == "WITH_KNOTS";
    KnotVector KU = knots(id, kind, pu, nu, wk ? kn[0] : nullptr, wk ? kn[2] : nullptr);
    KnotVector KV = knots(id, kind, pv, nv, wk ? kn[1] : nullptr, wk ? kn[3] : nullptr);
    std::vector<double> W;
    if (w)
      for (const auto& row : items(*w, id))
        for (const auto& x : items(row, id)) W.push_back(num(x, id));
    if (!W.empty() && W.size() != P.size()) bad(id, "number of weights differs from control points");
    return NurbsSurface(pu, pv, std::move(KU), std::move(KV), nu, nv, std::move(P), std::move(W)).clamped();
  }

  Curve3 curve(std::size_t id, int depth = 0) {
    if (depth > 8) bad(id, "recursive curve definition");
    const Instance& in = inst(id);
    Curve3 c;
    if (in.size() == 1) {
      const Record& r = in[0];
      if (r.name == "SURFACE_CURVE" || r.name == "SEAM_CURVE" || r.name == "INTERSECTION_CURVE")
        return curve(ref(arg(id, r, 1), id), depth + 1);
      if (r.name == "TRIMMED_CURVE") {  // the edge vertices bound the curve anyway
        c = curve(ref(arg(id, r, 1), id), depth + 1);
        if (!flag(arg(id, r, 4), id)) c.reversed = !c.reversed;
        return c;
      }
      if (r.name == "LINE") {
        c.kind = Curve3::Kind::Line;
        c.p = point(ref(arg(id, r, 1), id));
        const std::size_t vid = ref(arg(id, r, 2), id);
        const Record& vr = simple(vid, {"VECTOR"});
        c.d = direction(ref(arg(vid, vr, 1), vid));
        return c;
      }
      if (r.name == "CIRCLE" || r.name == "ELLIPSE") {
        c.kind = Curve3::Kind::Conic;
        c.f = placement(ref(arg(id, r, 1), id));
        c.a = len(arg(id, r, 2), id);
        c.b = r.name == "CIRCLE" ? c.a : len(arg(id, r, 3), id);
        if (!(c.a > 0.0) || !(c.b > 0.0)) bad(id, "non-positive radius");
        return c;
      }
      if (r.name == "HYPERBOLA" || r.name == "PARABOLA") {  // C + a cosh t X + b sinh t Y; C + a t^2 X + 2 a t Y
        c.kind = r.name == "HYPERBOLA" ? Curve3::Kind::Hyperbola : Curve3::Kind::Parabola;
        c.f = placement(ref(arg(id, r, 1), id));
        c.a = len(arg(id, r, 2), id);
        c.b = c.kind == Curve3::Kind::Hyperbola ? len(arg(id, r, 3), id) : c.a;
        if (!(c.a > 0.0) || !(c.b > 0.0)) bad(id, "non-positive conic parameter");
        return c;
      }
      if (r.name == "POLYLINE") {
        c.kind = Curve3::Kind::Spline;
        c.spline = make_polyline(points(arg(id, r, 1), id));
        return c;
      }
    }
    if (!bspline_kind(in, "B_SPLINE_CURVE").empty()) {
      c.kind = Curve3::Kind::Spline;
      c.spline = bspline_curve(id);
      return c;
    }
    unsupported(id);
  }

  Surf3 surface(std::size_t id) {
    const Instance& in = inst(id);
    Surf3 s;
    if (in.size() == 1) {
      const Record& r = in[0];
      const auto frame = [&] { return placement(ref(arg(id, r, 1), id)); };
      if (r.name == "PLANE") {
        s.kind = Surf3::Kind::Plane;
        s.f = frame();
        return s;
      }
      if (r.name == "CYLINDRICAL_SURFACE" || r.name == "SPHERICAL_SURFACE") {
        s.kind = r.name[0] == 'C' ? Surf3::Kind::Cylinder : Surf3::Kind::Sphere;
        s.f = frame();
        s.r1 = len(arg(id, r, 2), id);
        if (!(s.r1 > 0.0)) bad(id, "non-positive radius");
        return s;
      }
      if (r.name == "CONICAL_SURFACE") {
        s.kind = Surf3::Kind::Cone;
        s.f = frame();
        s.r1 = len(arg(id, r, 2), id);
        s.r2 = num(arg(id, r, 3), id) * angle_unit_;
        if (s.r1 < 0.0 || !(std::abs(s.r2) > 1e-12) || std::abs(s.r2) >= 0.5 * kPi) bad(id, "invalid cone");
        return s;
      }
      if (r.name == "TOROIDAL_SURFACE") {
        s.kind = Surf3::Kind::Torus;
        s.f = frame();
        s.r1 = len(arg(id, r, 2), id);
        s.r2 = len(arg(id, r, 3), id);
        if (!(s.r2 > 0.0) || !(s.r1 > s.r2)) bad(id, "unsupported torus radii");
        return s;
      }
    }
    if (!bspline_kind(in, "B_SPLINE_SURFACE").empty()) {
      s.kind = Surf3::Kind::Spline;
      s.spline = bspline_surface(id);
      return s;
    }
    unsupported(id);
  }

  // ---------------------------------------------------------------- edges

  /// The 3D curve of an EDGE_CURVE from its start to its end vertex.
  const NurbsCurve& edge(std::size_t id) {
    if (const auto it = edges_.find(id); it != edges_.end()) return it->second;
    const Record& r = simple(id, {"EDGE_CURVE"});
    const Vec3 ps = vertex(ref(arg(id, r, 1), id)), pe = vertex(ref(arg(id, r, 2), id));
    const Curve3 g = curve(ref(arg(id, r, 3), id));
    const bool same = flag(arg(id, r, 4), id) != g.reversed;
    const NurbsCurve c = with_end(with_start(edge_geometry(id, g, ps, pe, same), ps), pe);
    return edges_.emplace(id, c).first->second;
  }

  /// Edge curve from ps to pe on the curve g (same: edge direction = curve direction).
  NurbsCurve edge_geometry(std::size_t id, const Curve3& g, const Vec3& ps, const Vec3& pe, bool same) const {
    const bool closed = distance(ps, pe) <= tol_;
    switch (g.kind) {
      case Curve3::Kind::Line:
        if (closed) bad(id, "closed edge on a line");
        return make_line(ps, pe);
      case Curve3::Kind::Conic: {
        const auto ang = [&](const Vec3& p) {
          const Vec3 d = p - g.f.o;
          return std::atan2(dot(d, g.f.y) / g.b, dot(d, g.f.x) / g.a);
        };
        const double t0 = same ? ang(ps) : ang(pe), t1 = same ? ang(pe) : ang(ps);
        double sweep = closed ? kTwoPi : wrap_angle(t1 - t0);
        if (sweep <= 1e-12) sweep = kTwoPi;
        const NurbsCurve arc = make_arc({0, 0, 0}, {1, 0, 0}, {0, 1, 0}, 1.0, t0, t0 + sweep);
        const Vec3 X = g.a * g.f.x, Y = g.b * g.f.y, Z = g.f.z, O = g.f.o;
        const NurbsCurve c = arc.transformed(Transform({X.x, Y.x, Z.x, O.x, X.y, Y.y, Z.y, O.y, X.z, Y.z, Z.z, O.z}));
        return same ? c : c.reversed();
      }
      case Curve3::Kind::Hyperbola:
      case Curve3::Kind::Parabola: {
        if (closed) bad(id, "closed edge on an open conic");
        const bool hyp = g.kind == Curve3::Kind::Hyperbola;
        const auto param = [&](const Vec3& p) {
          const double y = dot(p - g.f.o, g.f.y);
          return hyp ? std::asinh(y / g.b) : y / (2.0 * g.a);
        };
        const auto at = [&](double t) {
          return hyp ? g.f.o + g.a * std::cosh(t) * g.f.x + g.b * std::sinh(t) * g.f.y
                     : g.f.o + g.a * t * t * g.f.x + 2.0 * g.a * t * g.f.y;
        };
        const auto deriv = [&](double t) {
          return hyp ? g.a * std::sinh(t) * g.f.x + g.b * std::cosh(t) * g.f.y : 2.0 * g.a * t * g.f.x + 2.0 * g.a * g.f.y;
        };
        double t0 = param(ps), t1 = param(pe);
        const bool rev = t0 > t1;  // edge against the curve direction (or inconsistent sense: follow the vertices)
        if (rev) std::swap(t0, t1);
        if (!(t1 - t0 > 1e-12)) bad(id, "degenerate edge");
        const NurbsCurve c = conic_arc(at, deriv, t0, t1, hyp);
        return rev ? c.reversed() : c;
      }
      case Curve3::Kind::Spline: break;
    }
    if (!g.spline) bad(id, "missing curve geometry");
    return spline_edge(*g.spline, ps, pe, same, closed, id);
  }

  NurbsCurve spline_edge(const NurbsCurve& c, const Vec3& ps, const Vec3& pe, bool same, bool closed_edge,
                         std::size_t id) const {
    const double a = c.domain_start(), b = c.domain_end(), eps = 1e-10 * (b - a);
    const bool closed_curve = distance(c.evaluate(a), c.evaluate(b)) <= tol_;
    double ta = closest_point(c, same ? ps : pe).u, tb = closest_point(c, same ? pe : ps).u;
    NurbsCurve out = c;
    if (closed_curve) {
      if (ta >= b - eps) ta = a;
      if (tb <= a + eps) tb = b;
      if (closed_edge) {
        if (ta > a + eps && ta < b - eps) out = join(c.sub_curve(ta, b), c.sub_curve(a, ta), tol_);
      } else if (tb > ta) {
        if (ta > a + eps || tb < b - eps) out = c.sub_curve(ta, tb);
      } else {
        out = join(c.sub_curve(ta, b), c.sub_curve(a, tb), tol_);
      }
    } else {
      if (closed_edge) bad(id, "closed edge on an open curve");
      const bool rev = ta > tb;
      if (rev) std::swap(ta, tb);
      if (tb - ta <= eps) bad(id, "degenerate edge");
      if (ta > a + eps || tb < b - eps) out = c.sub_curve(ta, tb);
      if (rev) out = out.reversed();
    }
    return same ? out : out.reversed();
  }

  /// Curve of an ORIENTED_EDGE in loop direction, and its EDGE_CURVE instance.
  std::pair<NurbsCurve, std::size_t> oriented_edge(std::size_t id, int depth = 0) {
    if (depth > 8) bad(id, "recursive oriented edge");
    const Record& r = simple(id, {"ORIENTED_EDGE", "EDGE_CURVE"});
    if (r.name == "EDGE_CURVE") return {edge(id), id};
    const bool fwd = flag(arg(id, r, 4), id);
    auto [c, e] = oriented_edge(ref(arg(id, r, 3), id), depth + 1);
    return {fwd ? c : c.reversed(), e};
  }

  // ---------------------------------------------------------------- faces

  std::vector<Loop3> face_loops(std::size_t id, const Record& face) {
    std::vector<Loop3> loops;
    for (const auto& bv : items(arg(id, face, 1), id)) {
      const std::size_t bid = ref(bv, id);
      const Record& br = simple(bid, {"FACE_BOUND", "FACE_OUTER_BOUND"});
      const std::size_t lid = ref(arg(bid, br, 1), bid);
      const bool fwd = flag(arg(bid, br, 2), bid);
      const Record& lr = simple(lid, {"EDGE_LOOP", "VERTEX_LOOP", "POLY_LOOP"});
      if (lr.name == "VERTEX_LOOP") continue;
      Loop3 l;
      if (lr.name == "POLY_LOOP") {
        const auto P = points(arg(lid, lr, 1), lid);
        for (std::size_t i = 0; i < P.size(); ++i) {
          l.curves.push_back(make_line(P[i], P[(i + 1) % P.size()]));
          l.edges.push_back(0);
        }
      } else {
        for (const auto& ov : items(arg(lid, lr, 1), lid)) {
          auto [c, e] = oriented_edge(ref(ov, lid));
          l.curves.push_back(std::move(c));
          l.edges.push_back(e);
        }
      }
      if (l.curves.empty()) continue;
      if (!fwd) {
        std::reverse(l.curves.begin(), l.curves.end());
        std::reverse(l.edges.begin(), l.edges.end());
        for (auto& c : l.curves) c = c.reversed();
      }
      loops.push_back(std::move(l));
    }
    return loops;
  }

  static Vec3 radial(const Frame& f, const Vec3& p) {
    const Vec3 d = p - f.o;
    const Vec3 r = d - dot(d, f.z) * f.z;
    return unit(r, f.x);
  }

  /// Normal of the elementary surface (S_u x S_v of its standard parametrisation).
  static Vec3 analytic_normal(const Surf3& s, const Vec3& p) {
    switch (s.kind) {
      case Surf3::Kind::Plane: return s.f.z;
      case Surf3::Kind::Cylinder: return radial(s.f, p);
      case Surf3::Kind::Sphere: return unit(p - s.f.o, s.f.z);
      case Surf3::Kind::Cone: return unit(radial(s.f, p) - std::tan(s.r2) * s.f.z, s.f.z);
      case Surf3::Kind::Torus: {
        const Vec3 c = s.f.o + s.r1 * radial(s.f, p);
        return unit(p - c, radial(s.f, p));
      }
      case Surf3::Kind::Spline: break;
    }
    return s.f.z;
  }

  /// Exact NURBS patch of an elementary surface covering the face boundary.
  NurbsSurface analytic_patch(std::size_t id, const Surf3& s, const std::vector<Loop3>& loops, bool same_sense) {
    struct BS {
      Vec3 p, inward;
      bool ok;
    };
    std::vector<std::vector<BS>> samples;  // one polyline per boundary curve
    std::vector<std::size_t> edge_of;
    std::map<std::size_t, int> uses;
    for (const auto& l : loops)
      for (std::size_t k = 0; k < l.curves.size(); ++k) {
        const NurbsCurve& c = l.curves[k];
        if (l.edges[k]) ++uses[l.edges[k]];
        edge_of.push_back(l.edges[k]);
        const double a = c.domain_start(), b = c.domain_end();
        const std::size_t per = std::max<std::size_t>(4, 64 / span_count(c, a, b));
        std::vector<BS> poly;
        for (double t : sample_params(c, a, b, per)) {
          const auto d = c.derivatives(t, 1);
          const Vec3 n = analytic_normal(s, d[0]) * (same_sense ? 1.0 : -1.0);
          const Vec3 in = cross(n, d[1]);
          const bool ok = norm(in) > 1e-9 * norm(d[1]) && norm(in) > 0.0;
          poly.push_back({d[0], ok ? in / norm(in) : Vec3{}, ok});
        }
        samples.push_back(std::move(poly));
      }
    const Frame& f = s.f;
    const auto theta = [&](const Vec3& p) {
      const Vec3 d = p - f.o;
      return std::atan2(dot(d, f.y), dot(d, f.x));
    };
    const auto e_theta = [&](const Vec3& p) { return cross(f.z, radial(f, p)); };
    const double size = std::max({s.r1, std::abs(s.r2), tol_});
    // angular coordinate samples: theta around the axis (phi around the tube of a torus)
    const auto ang_samples = [&](bool tube) {
      std::vector<std::vector<AngPt>> out;
      for (const auto& poly : samples) {
        std::vector<AngPt> a;
        for (const auto& q : poly) {
          const Vec3 d = q.p - f.o;
          const double rho = norm(d - dot(d, f.z) * f.z);
          AngPt ap;
          if (!tube) {
            ap.ok = rho > 1e-9 * size;
            ap.a = theta(q.p);
            ap.s = dot(q.inward, e_theta(q.p));
          } else {
            const double x = rho - s.r1, z = dot(d, f.z);
            ap.ok = std::hypot(x, z) > 1e-9 * size;
            ap.a = std::atan2(z, x);
            ap.s = dot(q.inward, -std::sin(ap.a) * radial(f, q.p) + std::cos(ap.a) * f.z);
          }
          ap.ok = ap.ok && q.ok;
          a.push_back(ap);
        }
        out.push_back(std::move(a));
      }
      return out;
    };
    // seam angle: a closed-surface seam edge (used twice by the face) of constant angle
    const auto range = [&](bool tube) {
      const auto a = ang_samples(tube);
      std::optional<double> seam;
      double fallback = 0.0;
      bool have_fallback = false;
      for (std::size_t i = 0; i < a.size(); ++i) {
        for (const auto& q : a[i])
          if (q.ok && !have_fallback) {
            fallback = q.a;
            have_fallback = true;
          }
        if (seam || !edge_of[i] || uses[edge_of[i]] < 2 || a[i].empty() || !a[i][0].ok) continue;
        const bool constant = std::all_of(a[i].begin(), a[i].end(), [&](const AngPt& q) {
          return !q.ok || std::abs(wrap_signed(q.a - a[i][0].a)) <= 1e-7;
        });
        if (constant) seam = a[i][0].a;
      }
      return angular_range(a, seam, fallback);
    };
    // axial range of a cylinder / cone
    const auto axial = [&] {
      double lo = 1e300, hi = -1e300;
      for (const auto& poly : samples)
        for (const auto& q : poly) {
          lo = std::min(lo, dot(q.p - f.o, f.z));
          hi = std::max(hi, dot(q.p - f.o, f.z));
        }
      if (!(hi - lo > tol_)) bad(id, "face boundary does not bound a region of the surface");
      const double m = 0.05 * (hi - lo);
      return std::make_pair(lo - m, hi + m);
    };
    const bool unbounded = samples.empty();
    const auto make_patch = [&]() -> NurbsSurface {
      switch (s.kind) {
        case Surf3::Kind::Plane: {
          if (unbounded) bad(id, "unbounded planar face");
          double x0 = 1e300, x1 = -1e300, y0 = 1e300, y1 = -1e300;
          for (const auto& poly : samples)
            for (const auto& q : poly) {
              const double x = dot(q.p - f.o, f.x), y = dot(q.p - f.o, f.y);
              x0 = std::min(x0, x);
              x1 = std::max(x1, x);
              y0 = std::min(y0, y);
              y1 = std::max(y1, y);
            }
          const double m = 0.05 * std::max(x1 - x0, y1 - y0);
          if (!(m > 0.0)) bad(id, "degenerate planar face");
          x0 -= m;
          x1 += m;
          y0 -= m;
          y1 += m;
          return make_rectangle(f.o + x0 * f.x + y0 * f.y, (x1 - x0) * f.x, (y1 - y0) * f.y);
        }
        case Surf3::Kind::Cylinder:
        case Surf3::Kind::Cone: {
          if (unbounded) bad(id, "unbounded face");
          const auto [t0, sweep] = range(false);
          auto [h0, h1] = axial();
          const Vec3 e = std::cos(t0) * f.x + std::sin(t0) * f.y;
          const double ta = s.kind == Surf3::Kind::Cone ? std::tan(s.r2) : 0.0;
          const auto rad = [&](double h) { return s.r1 + h * ta; };
          if (s.kind == Surf3::Kind::Cone) {  // do not extend beyond the apex
            const double apex = -s.r1 / ta;
            if (rad(h0) < 0.0) h0 = apex;
            if (rad(h1) < 0.0) h1 = apex;
          }
          const auto pt = [&](double h) {
            const double r = rad(h);
            return f.o + h * f.z + (std::abs(r) <= 1e-12 * size ? 0.0 : r) * e;
          };
          return revolve(make_line(pt(h0), pt(h1)), f.o, f.z, sweep);
        }
        case Surf3::Kind::Sphere: {
          const auto [t0, sweep] = unbounded ? std::make_pair(0.0, kTwoPi) : range(false);
          const Vec3 e = std::cos(t0) * f.x + std::sin(t0) * f.y;
          return revolve(make_arc(f.o, -f.z, e, s.r1, 0.0, kPi), f.o, f.z, sweep);
        }
        case Surf3::Kind::Torus: {
          const auto [t0, sweep] = unbounded ? std::make_pair(0.0, kTwoPi) : range(false);
          const auto [p0, psweep] = unbounded ? std::make_pair(0.0, kTwoPi) : range(true);
          const Vec3 e = std::cos(t0) * f.x + std::sin(t0) * f.y;
          return revolve(make_arc(f.o + s.r1 * e, e, f.z, s.r2, p0, p0 + psweep), f.o, f.z, sweep);
        }
        case Surf3::Kind::Spline: break;
      }
      bad(id, "not an elementary surface");
    };
    NurbsSurface out = make_patch();
    // orientation: S_u x S_v along the normal of the elementary surface
    const auto du = out.domain_u(), dv = out.domain_v();
    const double um = 0.5 * (du[0] + du[1]), vm = 0.5 * (dv[0] + dv[1]);
    if (dot(out.normal(um, vm), analytic_normal(s, out.evaluate(um, vm))) < 0.0) out = out.transposed();
    return out;
  }

  TrimmedSurface face(std::size_t id, bool shell_reversed, std::vector<std::string>& notes) {
    const Record& r = simple(id, {"ADVANCED_FACE", "FACE_SURFACE"});
    const bool same_sense = flag(arg(id, r, 3), id);
    const Surf3 s = surface(ref(arg(id, r, 2), id));
    const auto loops = face_loops(id, r);
    const NurbsSurface S = s.kind == Surf3::Kind::Spline ? *s.spline : analytic_patch(id, s, loops, same_sense);
    FaceMapper mapper(S, same_sense, tol_, notes);
    TrimmedSurface t = mapper.build(loops);
    return same_sense == shell_reversed ? t.flipped() : t;
  }

  void shell_faces(std::size_t id, bool reversed, std::vector<TrimmedSurface>& out) {
    const Record& r = simple(id, {"CLOSED_SHELL", "OPEN_SHELL", "ORIENTED_CLOSED_SHELL", "ORIENTED_OPEN_SHELL"});
    if (r.name.rfind("ORIENTED_", 0) == 0) {
      shell_faces(ref(arg(id, r, 2), id), reversed != !flag(arg(id, r, 3), id), out);
      return;
    }
    for (const auto& fv : items(arg(id, r, 1), id)) {
      const std::size_t fid = ref(fv, id);
      std::vector<std::string> notes;
      try {
        out.push_back(face(fid, reversed, notes));
      } catch (const std::exception& e) {
        problem(std::string(e.what()) + " (face #" + std::to_string(fid) + ")");
      }
      if (opt_.warnings)
        for (const auto& n : notes) opt_.warnings->push_back(n + " (face #" + std::to_string(fid) + ")");
    }
  }

  void add_shell(std::vector<Shell>& out, std::vector<TrimmedSurface> faces) const {
    if (faces.empty()) return;
    Shell sh = Shell::sew(std::move(faces), tol_);
    if (!sh.is_consistently_oriented()) sh = sh.oriented();
    out.push_back(std::move(sh));
  }
};

}  // namespace

std::vector<Shell> from_step(const std::string& text, const StepReadOptions& options) {
  try {
    Reader reader(Parser(text).parse(), options);
    return reader.read();
  } catch (const GeometryError& e) {
    const std::string what = e.what();
    if (what.rfind("read_step", 0) == 0) throw;
    throw GeometryError("read_step: " + what);
  }
}

std::vector<Shell> read_step(const std::string& path, const StepReadOptions& options) {
  std::ifstream f(path, std::ios::binary);
  if (!f) throw std::runtime_error("read_step: cannot open '" + path + "'");
  std::ostringstream ss;
  ss << f.rdbuf();
  return from_step(ss.str(), options);
}

}  // namespace nf
