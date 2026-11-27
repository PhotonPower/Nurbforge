#include "nurbforge/modeling/unroll.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <optional>
#include <string>
#include <vector>

#include "fit_family.hpp"
#include "nurbforge/analysis/properties.hpp"
#include "nurbforge/core/error.hpp"
#include "nurbforge/nurbs/join.hpp"

namespace nf {

namespace {

double size_of(const NurbsSurface& s) {
  double size = 0.0;
  for (const auto& p : s.control_points())
    size = std::max(size, distance(p, s.control_points().front()));
  return size;
}

// True if the surface is linear in v: S(u, v) lies on the line from S(u, v0) to S(u, v1) at the
// fraction of v.
bool linear_in_v(const NurbsSurface& s) {
  const auto du = s.domain_u(), dv = s.domain_v();
  const double eps = 1e-9 * std::max(size_of(s), 1e-300);
  const auto ku = unique_knots(s.knots_u()).first;
  for (std::size_t i = 0; i + 1 < ku.size(); ++i) {
    if (ku[i + 1] <= du[0] || ku[i] >= du[1]) continue;
    for (int q = 0; q <= 8; ++q) {
      const double u = ku[i] + (ku[i + 1] - ku[i]) * q / 8.0;
      const Vec3 a = s.evaluate(u, dv[0]), b = s.evaluate(u, dv[1]);
      for (const double f : {0.137, 0.5, 0.821})
        if (distance(s.evaluate(u, dv[0] + f * (dv[1] - dv[0])), a + f * (b - a)) > eps)
          return false;
    }
  }
  return true;
}

// Largest angle between the normals along the rulings u = const of a surface linear in v.
double twist_of(const NurbsSurface& s) {
  const auto du = s.domain_u(), dv = s.domain_v();
  const auto ku = unique_knots(s.knots_u()).first;
  double worst = 0.0;
  for (std::size_t i = 0; i + 1 < ku.size(); ++i) {
    if (ku[i + 1] <= du[0] || ku[i] >= du[1]) continue;
    for (int q = 0; q <= 16; ++q) {
      const double u = ku[i] + (ku[i + 1] - ku[i]) * q / 16.0;
      const Vec3 ref = s.normal(u, 0.5 * (dv[0] + dv[1]));
      for (int r = 0; r <= 8; ++r) {
        const Vec3 n = s.normal(u, dv[0] + (dv[1] - dv[0]) * r / 8.0);
        worst = std::max(worst, std::atan2(norm(cross(ref, n)), dot(ref, n)));
      }
    }
  }
  return worst;
}

constexpr int kGauss = 8;
constexpr double kNode[kGauss] = {-0.9602898564975363, -0.7966664774136267, -0.5255324099163290,
                                  -0.1834346424956498, 0.1834346424956498,  0.5255324099163290,
                                  0.7966664774136267,  0.9602898564975363};
constexpr double kWeight[kGauss] = {0.1012285362903763, 0.2223810344533745, 0.3137066458778873,
                                    0.3626837833783620, 0.3626837833783620, 0.3137066458778873,
                                    0.2223810344533745, 0.1012285362903763};

// Development of a surface linear in v: the directrix C(u) = S(u, vm) is laid out in the plane
// with its speed |C'| and its geodesic curvature, the rulings at their angles to it.
class Development {
 public:
  Development(const NurbsSurface& s, std::size_t sub) : s_(s) {
    const auto du = s.domain_u(), dv = s.domain_v();
    v0_ = dv[0];
    v1_ = dv[1];
    vm_ = 0.5 * (v0_ + v1_);
    const auto ku = unique_knots(s.knots_u()).first;
    grid_.push_back(du[0]);
    for (std::size_t i = 0; i + 1 < ku.size(); ++i) {
      if (ku[i + 1] <= du[0] || ku[i] >= du[1]) continue;
      for (std::size_t q = 1; q <= sub; ++q)
        grid_.push_back(q == sub ? ku[i + 1]
                                 : ku[i] + (ku[i + 1] - ku[i]) * static_cast<double>(q) /
                                               static_cast<double>(sub));
    }
    theta_.assign(grid_.size(), 0.0);
    pos_.assign(grid_.size(), Vec3{0, 0, 0});
    for (std::size_t i = 0; i + 1 < grid_.size(); ++i) {
      theta_[i + 1] = theta_[i] + turn(grid_[i], grid_[i + 1]);
      pos_[i + 1] = pos_[i] + advance(i, grid_[i + 1]);
    }
  }

  [[nodiscard]] const std::vector<double>& grid() const { return grid_; }
  [[nodiscard]] const std::vector<Vec3>& positions() const { return pos_; }

  // The images of S(u, v0) and S(u, v1).
  [[nodiscard]] std::array<Vec3, 2> rulings_ends(double u) const {
    const std::size_t i = interval(u);
    const double th = theta_[i] + turn(grid_[i], u);
    const Vec3 c = pos_[i] + advance(i, u);
    const Vec3 t{std::cos(th), std::sin(th), 0.0}, n{-std::sin(th), std::cos(th), 0.0};
    const Frame f = frame(u);
    const Vec3 D = s_.evaluate(u, v1_) - s_.evaluate(u, v0_);
    // the ruling in the tangent plane: components along the directrix and across it
    const Vec3 d = dot(D, f.t) * t + dot(D, cross(f.n, f.t)) * n;
    return {c - 0.5 * d, c + 0.5 * d};
  }

 private:
  struct Frame {
    Vec3 t, n;      // unit tangent of the directrix, surface normal
    double speed;   // |C'|
    double dtheta;  // geodesic curvature x speed: the turning rate of the developed directrix
  };

  [[nodiscard]] Frame frame(double u) const {
    const auto D = s_.derivatives(u, vm_, 2);
    const Vec3 c1 = D[1][0], c2 = D[2][0];
    const double speed = norm(c1);
    if (!(speed > 0.0))
      throw GeometryError(
          "unroll: the surface is degenerate along the middle of its rulings (u = " +
          std::to_string(u) + ")");
    const Vec3 n = s_.normal(u, vm_);
    return {c1 / speed, n, speed, dot(c2, cross(n, c1)) / (speed * speed)};
  }

  [[nodiscard]] std::size_t interval(double u) const {
    const auto it = std::upper_bound(grid_.begin(), grid_.end(), u);
    const auto i = static_cast<std::size_t>(std::max<std::ptrdiff_t>(0, it - grid_.begin() - 1));
    return std::min(i, grid_.size() - 2);
  }

  // change of the direction angle of the developed directrix from a to b (one Gauss rule)
  [[nodiscard]] double turn(double a, double b) const {
    double sum = 0.0;
    for (int k = 0; k < kGauss; ++k)
      sum += kWeight[k] * frame(0.5 * (a + b) + 0.5 * (b - a) * kNode[k]).dtheta;
    return 0.5 * (b - a) * sum;
  }

  // displacement of the developed directrix from grid point i to u
  [[nodiscard]] Vec3 advance(std::size_t i, double u) const {
    const double a = grid_[i];
    Vec3 sum{0, 0, 0};
    for (int k = 0; k < kGauss; ++k) {
      const double x = 0.5 * (a + u) + 0.5 * (u - a) * kNode[k];
      const double th = theta_[i] + turn(a, x);
      sum += kWeight[k] * frame(x).speed * Vec3{std::cos(th), std::sin(th), 0.0};
    }
    return 0.5 * (u - a) * sum;
  }

  const NurbsSurface& s_;
  double v0_{0.0}, v1_{1.0}, vm_{0.5};
  std::vector<double> grid_;
  std::vector<double> theta_;
  std::vector<Vec3> pos_;
};

// Pattern of a surface linear in v.
UnrolledSurface unroll_v(const NurbsSurface& s, double tolerance) {
  // quadrature: sub-intervals per knot span doubled until the developed directrix is settled
  std::optional<Development> dev;
  for (std::size_t sub = 2;; sub *= 2) {
    Development next(s, sub);
    if (dev) {
      double change = 0.0;
      for (std::size_t i = 0; i < dev->positions().size(); ++i)
        change = std::max(change, distance(dev->positions()[i], next.positions()[2 * i]));
      dev.emplace(std::move(next));
      if (change <= 0.05 * tolerance) break;
    } else {
      dev.emplace(std::move(next));
    }
    if (sub > 512)
      throw GeometryError("unroll: tolerance " + std::to_string(tolerance) +
                          " not reached by the quadrature");
  }
  const auto dv = s.domain_v();
  const auto exact = [&](double u) {
    const auto e = dev->rulings_ends(u);
    return std::vector<Vec3>{e[0], e[1]};
  };
  double max_error = 0.0;
  const auto error = [&](const std::vector<NurbsCurve>& c, double u, const std::vector<Vec3>& e) {
    const double err = std::max(distance(c[0].evaluate(u), e[0]), distance(c[1].evaluate(u), e[1]));
    if (err <= 0.9 * tolerance) max_error = std::max(max_error, err);
    return err;
  };
  // one piece per knot span of the surface (the development is smooth in between only)
  std::optional<NurbsSurface> pattern;
  const auto du = s.domain_u();
  std::vector<double> breaks{du[0]};
  for (const double k : unique_knots(s.knots_u()).first)
    if (k > breaks.back() + tol::kParametric && k <= du[1]) breaks.push_back(k);
  for (std::size_t i = 0; i + 1 < breaks.size(); ++i) {
    std::vector<double> T;
    for (int k = 0; k <= 4; ++k)
      T.push_back(k == 4 ? breaks[i + 1] : breaks[i] + (breaks[i + 1] - breaks[i]) * k / 4.0);
    const auto c = detail::fit_family(T, exact, error, 0.9 * tolerance, "unroll");
    std::vector<Vec3> P;
    for (std::size_t k = 0; k < c[0].num_control_points(); ++k) {
      P.push_back(c[0].control_points()[k]);
      P.push_back(c[1].control_points()[k]);
    }
    NurbsSurface piece{c[0].degree(),
                       1,
                       c[0].knots(),
                       KnotVector{dv[0], dv[0], dv[1], dv[1]},
                       c[0].num_control_points(),
                       2,
                       std::move(P)};
    pattern = pattern ? join_u(*pattern, piece, tolerance) : std::move(piece);
  }
  return {std::move(*pattern), max_error};
}

void measure(const NurbsSurface& s, bool ruled, Developability& d) {
  if (ruled) d.max_twist = twist_of(s);
  const auto du = s.domain_u(), dv = s.domain_v();
  for (int i = 0; i <= 12; ++i)
    for (int j = 0; j <= 12; ++j) {
      const double u = du[0] + (du[1] - du[0]) * (i + 0.5) / 13.0;
      const double v = dv[0] + (dv[1] - dv[0]) * (j + 0.5) / 13.0;
      try {
        d.max_gaussian_curvature =
            std::max(d.max_gaussian_curvature, std::abs(curvature(s, u, v).gaussian));
      } catch (const GeometryError&) {
        // no curvature at a degenerate point
      }
    }
}

}  // namespace

Developability developability(const NurbsSurface& surface, double angular_tolerance) {
  const NurbsSurface s = surface.clamped();
  Developability d;
  if (linear_in_v(s)) {
    d.rulings = RulingDirection::V;
    measure(s, true, d);
  } else if (const NurbsSurface t = s.transposed(); linear_in_v(t)) {
    d.rulings = RulingDirection::U;
    measure(t, true, d);
  } else {
    measure(s, false, d);
  }
  d.developable = d.rulings != RulingDirection::None && d.max_twist <= angular_tolerance;
  return d;
}

UnrolledSurface unroll(const NurbsSurface& surface, double tolerance, double angular_tolerance) {
  if (!(tolerance > 0.0)) throw GeometryError("unroll: tolerance must be positive");
  const Developability d = developability(surface, angular_tolerance);
  if (d.rulings == RulingDirection::None)
    throw GeometryError(
        "unroll: the surface is not ruled along a parameter direction (its parameter lines are "
        "not straight)");
  if (!d.developable)
    throw GeometryError("unroll: the surface is not developable (the normal turns by " +
                        std::to_string(d.max_twist) + " rad along a ruling)");
  const NurbsSurface s = surface.clamped();
  if (d.rulings == RulingDirection::V) return unroll_v(s, tolerance);
  // rulings along u: unroll the transposed surface; transposing back mirrors the orientation,
  // which the reflection y -> -y undoes
  UnrolledSurface r = unroll_v(s.transposed(), tolerance);
  std::vector<Vec3> P = r.pattern.control_points();
  for (auto& p : P) p.y = -p.y;
  const NurbsSurface mirrored{r.pattern.degree_u(), r.pattern.degree_v(), r.pattern.knots_u(),
                              r.pattern.knots_v(),  r.pattern.num_u(),    r.pattern.num_v(),
                              std::move(P)};
  return {mirrored.transposed(), r.max_error};
}

namespace {

struct Spring {
  std::size_t a, b;
  double length, weight;
};

// Geodesic distance of the neighbouring surface points at (ua, va) and (ub, vb): the chord,
// lengthened for the bulge of the surface between them. The bulge h is the distance of the
// surface from the chord's midpoint along the normal (an arc with sagitta h over the chord c
// has the length c (1 + 8 h^2 / (3 c^2))). Only the normal part counts: the parameter line
// between the points may bend within the surface, the geodesic does not.
double segment_length(const NurbsSurface& s, double ua, double va, double ub, double vb) {
  const Vec3 A = s.evaluate(ua, va), B = s.evaluate(ub, vb);
  const double c = distance(A, B);
  if (!(c > 0.0)) return 0.0;
  const double um = 0.5 * (ua + ub), vm = 0.5 * (va + vb);
  try {
    const double h = dot(s.evaluate(um, vm) - 0.5 * (A + B), s.normal(um, vm));
    return c * (1.0 + 8.0 * h * h / (3.0 * c * c));
  } catch (const GeometryError&) {
    return c;  // no normal at a degenerate point
  }
}

// Stress majorisation: positions X (x, y per node) minimising sum w (|x_a - x_b| - L)^2.
double relax(const std::vector<Spring>& springs, std::vector<std::array<double, 2>>& X,
             const FlattenOptions& opt) {
  const std::size_t n = X.size();
  std::vector<double> diag(n, 0.0);
  for (const auto& e : springs) {
    diag[e.a] += e.weight;
    diag[e.b] += e.weight;
  }
  using Field = std::vector<std::array<double, 2>>;
  const auto apply = [&](const Field& x, Field& y) {  // weighted graph Laplacian
    for (auto& v : y) v = {0.0, 0.0};
    for (const auto& e : springs)
      for (std::size_t k = 0; k < 2; ++k) {
        const double d = e.weight * (x[e.a][k] - x[e.b][k]);
        y[e.a][k] += d;
        y[e.b][k] -= d;
      }
  };
  const auto inner = [](const Field& a, const Field& b) {
    double sum = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) sum += a[i][0] * b[i][0] + a[i][1] * b[i][1];
    return sum;
  };
  const auto stress = [&](const Field& x) {
    double sum = 0.0;
    for (const auto& e : springs) {
      const double d = std::hypot(x[e.a][0] - x[e.b][0], x[e.a][1] - x[e.b][1]) - e.length;
      sum += e.weight * d * d;
    }
    return sum;
  };
  Field rhs(n), r(n), z(n), p(n), q(n);
  double last = stress(X);
  for (std::size_t it = 0; it < opt.max_iterations; ++it) {
    // right-hand side B(X) X
    for (auto& v : rhs) v = {0.0, 0.0};
    for (const auto& e : springs) {
      const double dx = X[e.a][0] - X[e.b][0], dy = X[e.a][1] - X[e.b][1];
      const double d = std::hypot(dx, dy);
      if (!(d > 1e-300)) continue;
      const double c = e.weight * e.length / d;
      rhs[e.a][0] += c * dx;
      rhs[e.a][1] += c * dy;
      rhs[e.b][0] -= c * dx;
      rhs[e.b][1] -= c * dy;
    }
    // L X = rhs by conjugate gradients with the diagonal as preconditioner, from the current X
    // (the constant null space is not excited: rhs and all residuals sum to zero)
    apply(X, q);
    for (std::size_t i = 0; i < n; ++i) r[i] = {rhs[i][0] - q[i][0], rhs[i][1] - q[i][1]};
    const double target = 1e-24 * std::max(inner(rhs, rhs), 1e-300);
    for (std::size_t i = 0; i < n; ++i) z[i] = {r[i][0] / diag[i], r[i][1] / diag[i]};
    p = z;
    double rz = inner(r, z);
    for (int cg = 0; cg < 400 && inner(r, r) > target; ++cg) {
      apply(p, q);
      const double alpha = rz / inner(p, q);
      for (std::size_t i = 0; i < n; ++i)
        for (std::size_t k = 0; k < 2; ++k) {
          X[i][k] += alpha * p[i][k];
          r[i][k] -= alpha * q[i][k];
        }
      for (std::size_t i = 0; i < n; ++i) z[i] = {r[i][0] / diag[i], r[i][1] / diag[i]};
      const double rz_new = inner(r, z);
      const double beta = rz_new / rz;
      rz = rz_new;
      for (std::size_t i = 0; i < n; ++i)
        for (std::size_t k = 0; k < 2; ++k) p[i][k] = z[i][k] + beta * p[i][k];
    }
    const double now = stress(X);
    const bool settled = last - now <= opt.convergence * std::max(last, 1e-300);
    last = now;
    // the majorisation is robust but converges linearly: a few steps, then Gauss-Newton
    if (settled || it + 1 >= std::min<std::size_t>(40, opt.max_iterations)) break;
  }

  // Gauss-Newton on the residuals sqrt(w) (|x_a - x_b| - L): quadratic convergence where the
  // surface is (nearly) developable. Steps that do not lower the stress are halved.
  Field grad(n), step(n), trial(n);
  std::vector<std::array<double, 2>> dir(springs.size());
  for (std::size_t it = 0; it < opt.max_iterations; ++it) {
    for (auto& v : grad) v = {0.0, 0.0};
    Field hd(n, {0.0, 0.0});
    for (std::size_t k = 0; k < springs.size(); ++k) {
      const Spring& e = springs[k];
      const double dx = X[e.a][0] - X[e.b][0], dy = X[e.a][1] - X[e.b][1];
      const double d = std::hypot(dx, dy);
      dir[k] = d > 1e-300 ? std::array<double, 2>{dx / d, dy / d} : std::array<double, 2>{0.0, 0.0};
      const double c = e.weight * (d - e.length);
      for (std::size_t q2 = 0; q2 < 2; ++q2) {
        grad[e.a][q2] += c * dir[k][q2];
        grad[e.b][q2] -= c * dir[k][q2];
        const double h = e.weight * dir[k][q2] * dir[k][q2];
        hd[e.a][q2] += h;
        hd[e.b][q2] += h;
      }
    }
    double hmean = 0.0;
    for (const auto& v : hd) hmean += v[0] + v[1];
    const double mu = 1e-9 * hmean / static_cast<double>(2 * n);  // lifts the rigid motions
    const auto hessian = [&](const Field& x, Field& y) {
      for (std::size_t i = 0; i < n; ++i) y[i] = {mu * x[i][0], mu * x[i][1]};
      for (std::size_t k = 0; k < springs.size(); ++k) {
        const Spring& e = springs[k];
        const double t =
            e.weight * (dir[k][0] * (x[e.a][0] - x[e.b][0]) + dir[k][1] * (x[e.a][1] - x[e.b][1]));
        for (std::size_t q2 = 0; q2 < 2; ++q2) {
          y[e.a][q2] += t * dir[k][q2];
          y[e.b][q2] -= t * dir[k][q2];
        }
      }
    };
    // H step = -grad by preconditioned conjugate gradients
    for (auto& v : step) v = {0.0, 0.0};
    for (std::size_t i = 0; i < n; ++i) r[i] = {-grad[i][0], -grad[i][1]};
    const double g2 = inner(grad, grad);
    if (!(g2 > 0.0)) break;
    const auto precondition = [&](const Field& x, Field& y) {
      for (std::size_t i = 0; i < n; ++i)
        y[i] = {x[i][0] / (hd[i][0] + mu), x[i][1] / (hd[i][1] + mu)};
    };
    precondition(r, z);
    p = z;
    double rz = inner(r, z);
    for (int cg = 0; cg < 600 && inner(r, r) > 1e-20 * g2; ++cg) {
      hessian(p, q);
      const double alpha = rz / inner(p, q);
      for (std::size_t i = 0; i < n; ++i)
        for (std::size_t k = 0; k < 2; ++k) {
          step[i][k] += alpha * p[i][k];
          r[i][k] -= alpha * q[i][k];
        }
      precondition(r, z);
      const double rz_new = inner(r, z);
      const double beta = rz_new / rz;
      rz = rz_new;
      for (std::size_t i = 0; i < n; ++i)
        for (std::size_t k = 0; k < 2; ++k) p[i][k] = z[i][k] + beta * p[i][k];
    }
    double t = 1.0, now = last;
    bool improved = false;
    for (int half = 0; half < 30; ++half, t *= 0.5) {
      for (std::size_t i = 0; i < n; ++i)
        trial[i] = {X[i][0] + t * step[i][0], X[i][1] + t * step[i][1]};
      now = stress(trial);
      if (now < last) {
        improved = true;
        break;
      }
    }
    if (!improved) break;
    X = trial;
    const bool settled = last - now <= opt.convergence * std::max(last, 1e-300) ||
                         now <= 1e-26 * static_cast<double>(springs.size());
    last = now;
    if (settled) break;
  }
  return last;
}

}  // namespace

FlattenedSurface flatten(const NurbsSurface& surface, const FlattenOptions& opt) {
  const std::size_t m = opt.samples_u, n = opt.samples_v;
  if (m < 4 || n < 4) throw GeometryError("flatten: at least 4 samples per direction are needed");
  const NurbsSurface s = surface.clamped();
  const auto du = s.domain_u(), dv = s.domain_v();
  std::vector<double> us(m), vs(n);
  for (std::size_t i = 0; i < m; ++i)
    us[i] = i + 1 == m
                ? du[1]
                : du[0] + (du[1] - du[0]) * static_cast<double>(i) / static_cast<double>(m - 1);
  for (std::size_t j = 0; j < n; ++j)
    vs[j] = j + 1 == n
                ? dv[1]
                : dv[0] + (dv[1] - dv[0]) * static_cast<double>(j) / static_cast<double>(n - 1);
  const auto id = [n](std::size_t i, std::size_t j) { return i * n + j; };

  // springs along the parameter lines and both diagonals of every cell
  std::vector<Spring> springs;
  std::vector<double> along_u((m - 1) * n), along_v(m * (n - 1)), across((m - 1) * (n - 1));
  double total = 0.0;
  for (std::size_t i = 0; i < m; ++i)
    for (std::size_t j = 0; j < n; ++j) {
      if (i + 1 < m) {
        const double L = segment_length(s, us[i], vs[j], us[i + 1], vs[j]);
        along_u[i * n + j] = L;
        springs.push_back({id(i, j), id(i + 1, j), L, 1.0});
      }
      if (j + 1 < n) {
        const double L = segment_length(s, us[i], vs[j], us[i], vs[j + 1]);
        along_v[i * (n - 1) + j] = L;
        springs.push_back({id(i, j), id(i, j + 1), L, 1.0});
      }
      if (i + 1 < m && j + 1 < n) {
        springs.push_back({id(i, j), id(i + 1, j + 1),
                           segment_length(s, us[i], vs[j], us[i + 1], vs[j + 1]), 1.0});
        const double L = segment_length(s, us[i + 1], vs[j], us[i], vs[j + 1]);
        across[i * (n - 1) + j] = L;
        springs.push_back({id(i + 1, j), id(i, j + 1), L, 1.0});
      }
    }
  for (const auto& e : springs) total += e.length;
  const double mean = total / static_cast<double>(springs.size());
  if (!(mean > 0.0)) throw GeometryError("flatten: the surface is degenerate");
  // relative errors: weight 1 / length (bounded for the collapsed edges at a pole)
  for (auto& e : springs) e.weight = 1.0 / std::max(e.length, 1e-3 * mean);

  // start: arc lengths along the parameter lines, measured from the middle of the grid
  const std::size_t ci = m / 2, cj = n / 2;
  std::vector<std::array<double, 2>> X(m * n);
  for (std::size_t j = 0; j < n; ++j) {
    double x = 0.0;
    X[id(ci, j)][0] = 0.0;
    for (std::size_t i = ci; i + 1 < m; ++i) X[id(i + 1, j)][0] = (x += along_u[i * n + j]);
    x = 0.0;
    for (std::size_t i = ci; i-- > 0;) X[id(i, j)][0] = (x -= along_u[i * n + j]);
  }
  for (std::size_t i = 0; i < m; ++i) {
    double y = 0.0;
    X[id(i, cj)][1] = 0.0;
    for (std::size_t j = cj; j + 1 < n; ++j) X[id(i, j + 1)][1] = (y += along_v[i * (n - 1) + j]);
    y = 0.0;
    for (std::size_t j = cj; j-- > 0;) X[id(i, j)][1] = (y -= along_v[i * (n - 1) + j]);
  }
  relax(springs, X, opt);

  // normal +z: the cells keep the orientation of (u, v); the middle at the origin, u along +x
  double signed_area = 0.0;
  for (std::size_t i = 0; i + 1 < m; ++i)
    for (std::size_t j = 0; j + 1 < n; ++j) {
      const auto& a = X[id(i, j)];
      const auto& b = X[id(i + 1, j + 1)];
      const auto& c = X[id(i + 1, j)];
      const auto& d = X[id(i, j + 1)];
      signed_area += 0.5 * ((b[0] - a[0]) * (d[1] - c[1]) - (b[1] - a[1]) * (d[0] - c[0]));
    }
  if (signed_area < 0.0)
    for (auto& x : X) x[1] = -x[1];
  const std::array<double, 2> origin = X[id(ci, cj)];
  const double ex = X[id(ci + 1, cj)][0] - X[id(ci - 1, cj)][0];
  const double ey = X[id(ci + 1, cj)][1] - X[id(ci - 1, cj)][1];
  const double el = std::hypot(ex, ey);
  const double co = el > 0.0 ? ex / el : 1.0, si = el > 0.0 ? ey / el : 0.0;
  for (auto& x : X) {
    const double px = x[0] - origin[0], py = x[1] - origin[1];
    x = {co * px + si * py, -si * px + co * py};
  }

  // pattern: bicubic interpolation of the grid at its parameters (rows along u, then along v)
  std::vector<std::vector<Vec3>> rows;  // rows[j]: control points along u for v = vs[j]
  KnotVector ku, kv;
  for (std::size_t j = 0; j < n; ++j) {
    std::vector<Vec3> pts(m);
    for (std::size_t i = 0; i < m; ++i) pts[i] = {X[id(i, j)][0], X[id(i, j)][1], 0.0};
    const NurbsCurve c = interpolate_curve(pts, 3, us);
    ku = c.knots();
    rows.push_back(c.control_points());
  }
  const std::size_t nu = rows.front().size();
  std::vector<Vec3> P;
  std::size_t nv = 0;
  for (std::size_t k = 0; k < nu; ++k) {
    std::vector<Vec3> pts(n);
    for (std::size_t j = 0; j < n; ++j) pts[j] = rows[j][k];
    const NurbsCurve c = interpolate_curve(pts, 3, vs);
    kv = c.knots();
    nv = c.control_points().size();
    P.insert(P.end(), c.control_points().begin(), c.control_points().end());
  }
  FlattenedSurface r{
      NurbsSurface{3, 3, ku, kv, nu, nv, std::move(P)}, {}, {}, {}, {}, 1.0, 1.0, 0.0, 0.0, 0.0};

  // Strains of the layout itself, cell by cell: each cell is two triangles (cut by the
  // diagonal from (i+1, j) to (i, j+1)) whose side lengths are known on the surface and in the
  // plane; the principal stretches are the singular values of the linear map between them.
  const auto stretches = [](double La, double Lb, double Lc, double la, double lb, double lc,
                            double& s1, double& s2) {
    // triangles with the corners (0, 0), (La, 0), (cx, cy) and likewise in the plane
    if (!(La > 0.0) || !(la > 0.0)) return false;
    const double cx = (La * La + Lb * Lb - Lc * Lc) / (2.0 * La);
    const double qx = (la * la + lb * lb - lc * lc) / (2.0 * la);
    const double cy2 = Lb * Lb - cx * cx, qy2 = lb * lb - qx * qx;
    if (!(cy2 > 1e-12 * Lb * Lb) || !(qy2 > 0.0)) return false;
    const double cy = std::sqrt(cy2), qy = std::sqrt(qy2);
    // F = [[p, q], [0, r]] maps the surface triangle onto the plane triangle
    const double p = la / La, q = (qx - p * cx) / cy, rr = qy / cy;
    const double tr = p * p + q * q + rr * rr, det = p * rr;
    const double disc = std::sqrt(std::max(0.0, 0.25 * tr * tr - det * det));
    s1 = std::sqrt(std::max(0.0, 0.5 * tr + disc));
    s2 = std::sqrt(std::max(0.0, 0.5 * tr - disc));
    return true;
  };
  const auto plane_distance = [&](std::size_t a, std::size_t b) {
    return std::hypot(X[a][0] - X[b][0], X[a][1] - X[b][1]);
  };
  for (std::size_t i = 0; i + 1 < m; ++i) r.u.push_back(0.5 * (us[i] + us[i + 1]));
  for (std::size_t j = 0; j + 1 < n; ++j) r.v.push_back(0.5 * (vs[j] + vs[j + 1]));
  double sum2 = 0.0;
  std::size_t count = 0;
  for (std::size_t i = 0; i + 1 < m; ++i)
    for (std::size_t j = 0; j + 1 < n; ++j) {
      const std::size_t p00 = id(i, j), p10 = id(i + 1, j), p01 = id(i, j + 1),
                        p11 = id(i + 1, j + 1);
      const double Lc = across[i * (n - 1) + j], lc = plane_distance(p10, p01);
      double hi = 1.0, lo = 1.0, s1 = 1.0, s2 = 1.0;
      bool any = false;
      if (stretches(along_u[i * n + j], along_v[i * (n - 1) + j], Lc, plane_distance(p00, p10),
                    plane_distance(p00, p01), lc, s1, s2)) {
        hi = s1;
        lo = s2;
        any = true;
      }
      if (stretches(along_u[i * n + j + 1], along_v[(i + 1) * (n - 1) + j], Lc,
                    plane_distance(p11, p01), plane_distance(p11, p10), lc, s1, s2)) {
        hi = any ? std::max(hi, s1) : s1;
        lo = any ? std::min(lo, s2) : s2;
        any = true;
      }
      if (any) {
        r.max_stretch = std::max(r.max_stretch, hi);
        r.min_stretch = std::min(r.min_stretch, lo);
        sum2 += (hi - 1.0) * (hi - 1.0) + (lo - 1.0) * (lo - 1.0);
        count += 2;
      }
      r.stretch_max.push_back(hi);
      r.stretch_min.push_back(lo);
    }
  r.rms_strain = count > 0 ? std::sqrt(sum2 / static_cast<double>(count)) : 0.0;
  r.area = area(s);
  r.pattern_area = area(r.pattern);
  return r;
}

}  // namespace nf
