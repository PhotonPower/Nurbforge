#include "nurbforge/modeling/sweep.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <string>
#include <vector>

#include "nurbforge/core/error.hpp"
#include "nurbforge/geometry/transform.hpp"
#include "nurbforge/modeling/construct.hpp"
#include "nurbforge/modeling/loft.hpp"
#include "nurbforge/nurbs/join.hpp"

namespace nf {

namespace {

// 3x3 matrix stored by columns.
struct Mat3 {
  Vec3 c0{1, 0, 0}, c1{0, 1, 0}, c2{0, 0, 1};
};

Vec3 mul(const Mat3& M, const Vec3& x) { return x.x * M.c0 + x.y * M.c1 + x.z * M.c2; }

// A * B^T for orthonormal frames (columns of A and B are the frame axes).
Mat3 relative(const Mat3& A, const Mat3& B) {
  // (A B^T) e_k = A (B^T e_k) = A * (row k of B) = A * (c0.k, c1.k, c2.k)
  const Vec3 r0{B.c0.x, B.c1.x, B.c2.x}, r1{B.c0.y, B.c1.y, B.c2.y}, r2{B.c0.z, B.c1.z, B.c2.z};
  return {mul(A, r0), mul(A, r1), mul(A, r2)};
}

Transform linear(const Mat3& R) {
  return Transform({R.c0.x, R.c1.x, R.c2.x, 0.0, R.c0.y, R.c1.y, R.c2.y, 0.0, R.c0.z, R.c1.z,
                    R.c2.z, 0.0});
}

Vec3 unit_tangent(const NurbsCurve& c, double v) {
  const Vec3 d = c.derivatives(v, 1)[1];
  const double n = norm(d);
  if (n <= tol::kLength) throw GeometryError("sweep: trajectory tangent vanishes at v = " + std::to_string(v));
  return d / n;
}

Vec3 perpendicular(const Vec3& t) {
  const Vec3 helper = std::abs(t.x) < 0.9 ? Vec3{1, 0, 0} : Vec3{0, 1, 0};
  return normalized(cross(t, helper));
}

// Rotation-minimising frame by the double reflection method (Wang, Jüttler, Zheng, Liu 2008),
// tabulated at N+1 uniform parameters; frames in between by one extra reflection step.
class RotationMinimizingFrame {
 public:
  explicit RotationMinimizingFrame(const NurbsCurve& c, std::size_t N = 4096)
      : c_(c), a_(c.domain_start()), h_((c.domain_end() - c.domain_start()) / static_cast<double>(N)) {
    x_.resize(N + 1);
    t_.resize(N + 1);
    r_.resize(N + 1);
    for (std::size_t i = 0; i <= N; ++i) {
      const double v = a_ + h_ * static_cast<double>(i);
      x_[i] = c.evaluate(v);
      t_[i] = unit_tangent(c, v);
    }
    r_[0] = perpendicular(t_[0]);
    for (std::size_t i = 0; i < N; ++i) r_[i + 1] = step(x_[i], t_[i], r_[i], x_[i + 1], t_[i + 1]);
  }

  [[nodiscard]] Mat3 frame(double v) const {
    const std::size_t N = x_.size() - 1;
    const double s = (v - a_) / h_;
    const auto i = static_cast<std::size_t>(std::clamp(std::floor(s), 0.0, static_cast<double>(N)));
    const Vec3 t = unit_tangent(c_, v);
    const Vec3 r = step(x_[i], t_[i], r_[i], c_.evaluate(v), t);
    return {t, r, cross(t, r)};
  }

 private:
  static Vec3 step(const Vec3& x0, const Vec3& t0, const Vec3& r0, const Vec3& x1, const Vec3& t1) {
    const Vec3 v1 = x1 - x0;
    const double c1 = dot(v1, v1);
    Vec3 rL = r0, tL = t0;
    if (c1 > 0.0) {
      rL = r0 - (2.0 / c1) * dot(v1, r0) * v1;
      tL = t0 - (2.0 / c1) * dot(v1, t0) * v1;
    }
    const Vec3 v2 = t1 - tL;
    const double c2 = dot(v2, v2);
    Vec3 r1 = c2 > 0.0 ? rL - (2.0 / c2) * dot(v2, rL) * v2 : rL;
    r1 = r1 - dot(r1, t1) * t1;  // remove round-off drift
    return normalized(r1);
  }

  const NurbsCurve& c_;
  double a_, h_;
  std::vector<Vec3> x_, t_, r_;
};

using Placement = std::function<Transform(double)>;

// Skins one smooth piece [a, b] through exactly placed copies of the profile (on [0, 1]),
// doubling the number of sections until the deviation from the ideal sweep at the midpoints is
// below the tolerance.
NurbsSurface skin_piece(const NurbsCurve& profile, const Placement& place, double a, double b,
                        const SweepOptions& opt) {
  std::size_t K = 5;
  while (true) {
    std::vector<double> v(K);
    std::vector<NurbsCurve> sections;
    for (std::size_t k = 0; k < K; ++k) {
      v[k] = k + 1 == K ? b : a + (b - a) * static_cast<double>(k) / static_cast<double>(K - 1);
      sections.push_back(profile.transformed(place(v[k])));
    }
    NurbsSurface s = loft(sections, 3, v);
    double err = 0.0;
    for (std::size_t k = 0; k + 1 < K; ++k) {
      const double vm = 0.5 * (v[k] + v[k + 1]);
      const Transform T = place(vm);
      for (int i = 0; i <= 16; ++i) {
        const double u = i / 16.0;
        err = std::max(err, distance(s.evaluate(u, vm), T.apply(profile.evaluate(u))));
      }
    }
    if (err <= opt.tolerance) return s;
    if (2 * K - 1 > opt.max_sections)
      throw GeometryError("sweep: tolerance " + std::to_string(opt.tolerance) + " not reached with " +
                          std::to_string(K) + " sections per piece (deviation " + std::to_string(err) + ")");
    K = 2 * K - 1;
  }
}

// The placement is only as smooth as the guide curves, which are typically C1/C2 at their knots
// (e.g. the rational circle). A cubic skin across such a kink converges like h^2 only, so the
// sweep is skinned piecewise between the distinct knots and the pieces are joined (they share
// the exactly placed boundary section).
NurbsSurface skin_adaptive(const NurbsCurve& profile, const Placement& place,
                           std::vector<double> breaks, const SweepOptions& opt) {
  std::sort(breaks.begin(), breaks.end());
  breaks.erase(std::unique(breaks.begin(), breaks.end(),
                           [](double x, double y) { return y - x <= tol::kParametric; }),
               breaks.end());
  NurbsSurface result = skin_piece(profile, place, breaks[0], breaks[1], opt);
  for (std::size_t i = 1; i + 1 < breaks.size(); ++i)
    result = join_v(result, skin_piece(profile, place, breaks[i], breaks[i + 1], opt), 1e-9);
  return result;
}

// Distinct knots of the curves inside [a, b], including a and b.
std::vector<double> break_points(const std::vector<const NurbsCurve*>& curves, double a, double b) {
  std::vector<double> x{a, b};
  for (const auto* c : curves)
    for (double k : unique_knots(c->knots()).first)
      if (k > a && k < b) x.push_back(k);
  return x;
}

}  // namespace

NurbsSurface sweep_translational(const NurbsCurve& profile_in, const NurbsCurve& trajectory) {
  const NurbsCurve profile = reparametrize(profile_in.clamped(), 0.0, 1.0);
  const Vec3 start = trajectory.evaluate(trajectory.domain_start());
  const std::size_t nu = profile.num_control_points(), nv = trajectory.num_control_points();
  std::vector<Vec3> P;
  std::vector<double> W;
  for (std::size_t i = 0; i < nu; ++i)
    for (std::size_t j = 0; j < nv; ++j) {
      P.push_back(profile.control_points()[i] + trajectory.control_points()[j] - start);
      W.push_back(profile.weights()[i] * trajectory.weights()[j]);
    }
  return {profile.degree(), trajectory.degree(), profile.knots(), trajectory.knots(), nu, nv,
          std::move(P), std::move(W)};
}

NurbsSurface sweep(const NurbsCurve& profile_in, const NurbsCurve& trajectory, const SweepOptions& opt) {
  if (opt.frame == SweepFrame::Fixed && opt.twist == 0.0) return sweep_translational(profile_in, trajectory);
  const NurbsCurve profile = reparametrize(profile_in.clamped(), 0.0, 1.0);
  const double a = trajectory.domain_start(), b = trajectory.domain_end();
  const Vec3 x0 = trajectory.evaluate(a);
  const RotationMinimizingFrame rmf(trajectory);
  const Mat3 F0 = rmf.frame(a);
  const Placement place = [&](double v) {
    Mat3 R;  // identity for SweepFrame::Fixed
    if (opt.frame == SweepFrame::RotationMinimizing) R = relative(rmf.frame(v), F0);
    Transform L = linear(R);
    if (opt.twist != 0.0)
      L = Transform::rotation({0, 0, 0}, unit_tangent(trajectory, v), opt.twist * (v - a) / (b - a)) * L;
    return Transform::translation(trajectory.evaluate(v)) * L * Transform::translation(-x0);
  };
  return skin_adaptive(profile, place, break_points({&trajectory}, a, b), opt);
}

NurbsSurface sweep2(const NurbsCurve& profile_in, const NurbsCurve& rail1_in, const NurbsCurve& rail2_in,
                    const SweepOptions& opt) {
  const NurbsCurve profile = reparametrize(profile_in.clamped(), 0.0, 1.0);
  const NurbsCurve rail1 = rail1_in.clamped();
  const double a = rail1.domain_start(), b = rail1.domain_end();
  const NurbsCurve rail2 = reparametrize(rail2_in.clamped(), a, b);
  const auto frame = [&](double v, double& scale) {
    const Vec3 d = rail2.evaluate(v) - rail1.evaluate(v);
    scale = norm(d);
    if (scale <= tol::kLength) throw GeometryError("sweep2: the rails touch at v = " + std::to_string(v));
    const Vec3 x = d / scale;
    const Vec3 t = rail1.derivatives(v, 1)[1];
    const Vec3 n = t - dot(t, x) * x;
    if (norm(n) <= tol::kAngular * norm(t))
      throw GeometryError("sweep2: rail1 tangent is parallel to the rail-to-rail direction");
    const Vec3 z = n / norm(n);
    return Mat3{x, cross(z, x), z};
  };
  double s0 = 0.0;
  const Mat3 M0 = frame(a, s0);
  const Vec3 x0 = rail1.evaluate(a);
  const Placement place = [&](double v) {
    double s = 0.0;
    const Mat3 M = frame(v, s);
    Mat3 R = relative(M, M0);
    const double f = s / s0;
    R = {f * R.c0, f * R.c1, f * R.c2};
    return Transform::translation(rail1.evaluate(v)) * linear(R) * Transform::translation(-x0);
  };
  return skin_adaptive(profile, place, break_points({&rail1, &rail2}, a, b), opt);
}

}  // namespace nf
