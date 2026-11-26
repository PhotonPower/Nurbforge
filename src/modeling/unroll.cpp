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

}  // namespace nf
