#include "nurbforge/modeling/sweep.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <numbers>
#include <string>
#include <vector>

#include "nurbforge/analysis/closest_point.hpp"
#include "nurbforge/core/error.hpp"
#include "nurbforge/geometry/primitives.hpp"
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
  return Transform(
      {R.c0.x, R.c1.x, R.c2.x, 0.0, R.c0.y, R.c1.y, R.c2.y, 0.0, R.c0.z, R.c1.z, R.c2.z, 0.0});
}

Vec3 unit_tangent(const NurbsCurve& c, double v) {
  const Vec3 d = c.derivatives(v, 1)[1];
  const double n = norm(d);
  if (n <= tol::kLength)
    throw GeometryError("sweep: trajectory tangent vanishes at v = " + std::to_string(v));
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
      : c_(c),
        a_(c.domain_start()),
        h_((c.domain_end() - c.domain_start()) / static_cast<double>(N)) {
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
using SectionAt = std::function<NurbsCurve(double)>;  // the exact section (on [0, 1]) at v

// Skins one smooth piece [a, b] through exact sections, doubling the number of sections until
// the deviation from the ideal sweep at the midpoints is below the tolerance.
NurbsSurface skin_piece(const SectionAt& section, double a, double b, const SweepOptions& opt) {
  std::size_t K = 5;
  while (true) {
    std::vector<double> v(K);
    std::vector<NurbsCurve> sections;
    for (std::size_t k = 0; k < K; ++k) {
      v[k] = k + 1 == K ? b : a + (b - a) * static_cast<double>(k) / static_cast<double>(K - 1);
      sections.push_back(section(v[k]));
    }
    NurbsSurface s = loft(sections, 3, v);
    double err = 0.0;
    for (std::size_t k = 0; k + 1 < K; ++k) {
      const double vm = 0.5 * (v[k] + v[k + 1]);
      const NurbsCurve ideal = section(vm);
      for (int i = 0; i <= 16; ++i) {
        const double u = i / 16.0;
        err = std::max(err, distance(s.evaluate(u, vm), ideal.evaluate(u)));
      }
    }
    if (err <= opt.tolerance) return s;
    if (2 * K - 1 > opt.max_sections)
      throw GeometryError("sweep: tolerance " + std::to_string(opt.tolerance) +
                          " not reached with " + std::to_string(K) +
                          " sections per piece (deviation " + std::to_string(err) + ")");
    K = 2 * K - 1;
  }
}

// The placement is only as smooth as the guide curves, which are typically C1/C2 at their knots
// (e.g. the rational circle). A cubic skin across such a kink converges like h^2 only, so the
// sweep is skinned piecewise between the distinct knots and the pieces are joined (they share
// the exactly placed boundary section).
NurbsSurface skin_adaptive(const SectionAt& section, std::vector<double> breaks,
                           const SweepOptions& opt) {
  std::sort(breaks.begin(), breaks.end());
  breaks.erase(std::unique(breaks.begin(), breaks.end(),
                           [](double x, double y) { return y - x <= tol::kParametric; }),
               breaks.end());
  NurbsSurface result = skin_piece(section, breaks[0], breaks[1], opt);
  for (std::size_t i = 1; i + 1 < breaks.size(); ++i)
    result = join_v(result, skin_piece(section, breaks[i], breaks[i + 1], opt), 1e-9);
  return result;
}

NurbsSurface skin_adaptive(const NurbsCurve& profile, const Placement& place,
                           std::vector<double> breaks, const SweepOptions& opt) {
  return skin_adaptive([&](double v) { return profile.transformed(place(v)); }, std::move(breaks),
                       opt);
}

// Inverse of an affine map (3x3 part inverted by its adjugate).
Transform inverse(const Transform& T) {
  const auto& m = T.matrix();
  const Vec3 c0{m[0], m[4], m[8]}, c1{m[1], m[5], m[9]}, c2{m[2], m[6], m[10]};
  const double det = dot(c0, cross(c1, c2));
  if (std::abs(det) <= 1e-300) throw GeometryError("sweep: singular placement");
  // rows of the inverse: (c1 x c2, c2 x c0, c0 x c1) / det
  const Vec3 r0 = cross(c1, c2) / det, r1 = cross(c2, c0) / det, r2 = cross(c0, c1) / det;
  const Vec3 t{m[3], m[7], m[11]};
  return Transform({r0.x, r0.y, r0.z, -dot(r0, t), r1.x, r1.y, r1.z, -dot(r1, t), r2.x, r2.y, r2.z,
                    -dot(r2, t)});
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
          std::move(P),     std::move(W)};
}

namespace {

// Placement of the profile (given at the trajectory start) at trajectory parameter v.
Placement sweep_placement(const NurbsCurve& trajectory, const RotationMinimizingFrame& rmf,
                          const SweepOptions& opt) {
  const double a = trajectory.domain_start(), b = trajectory.domain_end();
  const Vec3 x0 = trajectory.evaluate(a);
  const Mat3 F0 = rmf.frame(a);
  return [&trajectory, &rmf, opt, a, b, x0, F0](double v) {
    Mat3 R;  // identity for SweepFrame::Fixed
    if (opt.frame == SweepFrame::RotationMinimizing) R = relative(rmf.frame(v), F0);
    Transform L = linear(R);
    if (opt.twist != 0.0)
      L = Transform::rotation({0, 0, 0}, unit_tangent(trajectory, v),
                              opt.twist * (v - a) / (b - a)) *
          L;
    return Transform::translation(trajectory.evaluate(v)) * L * Transform::translation(-x0);
  };
}

}  // namespace

NurbsSurface sweep(const NurbsCurve& profile_in, const NurbsCurve& trajectory,
                   const SweepOptions& opt) {
  if (opt.frame == SweepFrame::Fixed && opt.twist == 0.0)
    return sweep_translational(profile_in, trajectory);
  const NurbsCurve profile = reparametrize(profile_in.clamped(), 0.0, 1.0);
  const double a = trajectory.domain_start(), b = trajectory.domain_end();
  const RotationMinimizingFrame rmf(trajectory);
  const Placement place = sweep_placement(trajectory, rmf, opt);
  return skin_adaptive(profile, place, break_points({&trajectory}, a, b), opt);
}

namespace {

// The profiles on [0, 1] and their parameters on the guide curve (where each profile comes
// closest to it, if not given), checked; returns the order of increasing parameter.
std::vector<std::size_t> arrange_profiles(const std::string& who,
                                          const std::vector<NurbsCurve>& profiles_in,
                                          const NurbsCurve& guide, std::vector<double>& params,
                                          std::vector<NurbsCurve>& profiles) {
  if (profiles_in.empty()) throw GeometryError(who + ": no profile");
  if (!params.empty() && params.size() != profiles_in.size())
    throw GeometryError(who + ": one parameter per profile is needed");
  const double a = guide.domain_start(), b = guide.domain_end();
  for (const auto& p : profiles_in) profiles.push_back(reparametrize(p.clamped(), 0.0, 1.0));
  if (params.empty())
    for (const auto& p : profiles) {
      double best = std::numeric_limits<double>::infinity(), t = a;
      for (int i = 0; i <= 64; ++i) {
        const CurvePoint cp = closest_point(guide, p.evaluate(i / 64.0));
        if (cp.distance < best) {
          best = cp.distance;
          t = cp.u;
        }
      }
      params.push_back(t);
    }
  std::vector<std::size_t> order(profiles.size());
  for (std::size_t i = 0; i < order.size(); ++i) order[i] = i;
  std::sort(order.begin(), order.end(),
            [&](std::size_t x, std::size_t y) { return params[x] < params[y]; });
  for (std::size_t i = 0; i < order.size(); ++i) {
    const double t = params[order[i]];
    if (t < a - tol::kParametric || t > b + tol::kParametric)
      throw GeometryError(who + ": a profile parameter lies outside the domain of the guide curve");
    if (i > 0 && !(t > params[order[i - 1]] + tol::kParametric))
      throw GeometryError(who + ": two profiles at the same parameter");
  }
  return order;
}

// Sweep of several profiles: each is moved back to the start by the inverse of its placement,
// the moved profiles are blended by a loft over v, and the blended shape is placed at v.
NurbsSurface skin_blended(const std::vector<NurbsCurve>& profiles,
                          const std::vector<std::size_t>& order, const std::vector<double>& params,
                          const Placement& place, std::vector<double> breaks, double a, double b,
                          const SweepOptions& opt) {
  if (profiles.size() == 1)
    return skin_adaptive(profiles[0].transformed(inverse(place(params[0]))), place,
                         std::move(breaks), opt);
  std::vector<NurbsCurve> local;
  std::vector<double> v;
  for (const std::size_t i : order) {
    local.push_back(profiles[i].transformed(inverse(place(params[i]))));
    v.push_back(params[i]);
  }
  const NurbsSurface shape = loft(local, std::min<std::size_t>(3, local.size() - 1), v);
  const SectionAt section = [&](double t) {
    return reparametrize(shape.iso_curve_v(std::clamp(t, v.front(), v.back())), 0.0, 1.0)
        .transformed(place(t));
  };
  for (const double t : unique_knots(shape.knots_v()).first)
    if (t > a && t < b) breaks.push_back(t);
  return skin_adaptive(section, std::move(breaks), opt);
}

}  // namespace

NurbsSurface sweep(const std::vector<NurbsCurve>& profiles_in, const NurbsCurve& trajectory,
                   std::vector<double> params, const SweepOptions& opt) {
  std::vector<NurbsCurve> profiles;
  const auto order = arrange_profiles("sweep", profiles_in, trajectory, params, profiles);
  const double a = trajectory.domain_start(), b = trajectory.domain_end();
  const RotationMinimizingFrame rmf(trajectory);
  const Placement place = sweep_placement(trajectory, rmf, opt);
  if (profiles.size() == 1)  // a single profile placed anywhere: move it back to the start
    return sweep(profiles[0].transformed(inverse(place(params[0]))), trajectory, opt);
  return skin_blended(profiles, order, params, place, break_points({&trajectory}, a, b), a, b, opt);
}

NurbsSurface rail_revolve(const NurbsCurve& profile_in, const NurbsCurve& rail,
                          const Vec3& axis_point, const Vec3& axis_dir, bool scale_height,
                          const SweepOptions& opt) {
  const double na = norm(axis_dir);
  if (na <= tol::kLength) throw GeometryError("rail_revolve: axis direction must be non-zero");
  const Vec3 ax = axis_dir / na;
  const NurbsCurve profile = reparametrize(profile_in.clamped(), 0.0, 1.0);
  const double a = rail.domain_start(), b = rail.domain_end();
  // the profile end on the rail (the one nearer to it)
  const Vec3 p0 = profile.evaluate(0.0), p1 = profile.evaluate(1.0);
  const Vec3 E = closest_point(rail, p0).distance <= closest_point(rail, p1).distance ? p0 : p1;
  const auto height = [&](const Vec3& p) { return dot(p - axis_point, ax); };
  const auto radial = [&](const Vec3& p) { return p - axis_point - height(p) * ax; };
  const Vec3 rE = radial(E);
  const double rhoE = norm(rE), hE = height(E);
  if (rhoE <= tol::kLength)
    throw GeometryError("rail_revolve: the profile end on the rail lies on the axis");
  if (scale_height && std::abs(hE) <= tol::kLength)
    throw GeometryError(
        "rail_revolve: scale_height needs a profile end off the axis point's "
        "height");
  const Vec3 e1 = rE / rhoE, e2 = cross(ax, e1);
  const auto raw_angle = [&](double t) {
    const Vec3 r = radial(rail.evaluate(t));
    return std::atan2(dot(r, e2), dot(r, e1));
  };
  // unwrapped angle table (continuous along the rail)
  constexpr std::size_t N = 2048;
  std::vector<double> table(N + 1);
  table[0] = raw_angle(a);
  for (std::size_t i = 1; i <= N; ++i) {
    const double t = a + (b - a) * static_cast<double>(i) / N;
    double d = raw_angle(t) - table[i - 1];
    d -= 2 * std::numbers::pi * std::round(d / (2 * std::numbers::pi));
    table[i] = table[i - 1] + d;
  }
  const auto angle = [&](double t) {
    const double s = (t - a) / (b - a) * N;
    const auto i = static_cast<std::size_t>(std::clamp(std::round(s), 0.0, static_cast<double>(N)));
    double d = raw_angle(t) - table[i];
    d -= 2 * std::numbers::pi * std::round(d / (2 * std::numbers::pi));
    return table[i] + d;
  };
  const SectionAt section = [&](double t) {
    const Vec3 R = rail.evaluate(t);
    const double rho = norm(radial(R));
    if (rho <= tol::kLength) throw GeometryError("rail_revolve: the rail touches the axis");
    const double k = rho / rhoE, kh = scale_height ? height(R) / hE : 1.0;
    // p -> O + k (radial part) + kh (axial part): linear part k (I - a a^T) + kh a a^T
    std::array<double, 12> m{};
    const double A[3] = {ax.x, ax.y, ax.z};
    const double O[3] = {axis_point.x, axis_point.y, axis_point.z};
    for (int r = 0; r < 3; ++r) {
      double tr = O[r];
      for (int c = 0; c < 3; ++c) {
        const double L = (r == c ? k : 0.0) + (kh - k) * A[r] * A[c];
        m[static_cast<std::size_t>(4 * r + c)] = L;
        tr -= L * O[c];
      }
      m[static_cast<std::size_t>(4 * r + 3)] = tr;
    }
    return profile.transformed(Transform::rotation(axis_point, ax, angle(t)) * Transform(m));
  };
  return skin_adaptive(section, break_points({&rail}, a, b), opt);
}

NurbsSurface sweep2(const NurbsCurve& profile, const NurbsCurve& rail1, const NurbsCurve& rail2,
                    const SweepOptions& opt) {
  return sweep2(std::vector<NurbsCurve>{profile}, rail1, rail2, {rail1.domain_start()},
                Sweep2Scaling::Uniform, opt);
}

NurbsSurface sweep2(const std::vector<NurbsCurve>& profiles_in, const NurbsCurve& rail1_in,
                    const NurbsCurve& rail2_in, std::vector<double> params, Sweep2Scaling scaling,
                    const SweepOptions& opt) {
  const NurbsCurve rail1 = rail1_in.clamped();
  const double a = rail1.domain_start(), b = rail1.domain_end();
  const NurbsCurve rail2 = reparametrize(rail2_in.clamped(), a, b);
  std::vector<NurbsCurve> profiles;
  const auto order = arrange_profiles("sweep2", profiles_in, rail1, params, profiles);
  const auto frame = [&](double v, double& scale) {
    const Vec3 d = rail2.evaluate(v) - rail1.evaluate(v);
    scale = norm(d);
    if (scale <= tol::kLength)
      throw GeometryError("sweep2: the rails touch at v = " + std::to_string(v));
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
    // M diag(f, g, g) M0^T: f along the rail-to-rail axis, g normal to it
    const double f = s / s0, g = scaling == Sweep2Scaling::Uniform ? f : 1.0;
    const Mat3 R = relative(Mat3{f * M.c0, g * M.c1, g * M.c2}, M0);
    return Transform::translation(rail1.evaluate(v)) * linear(R) * Transform::translation(-x0);
  };
  std::vector<double> breaks = break_points({&rail1, &rail2}, a, b);
  if (profiles.size() == 1 && params[0] == a)  // given at the rail starts: no moving back
    return skin_adaptive(profiles[0], place, std::move(breaks), opt);
  return skin_blended(profiles, order, params, place, std::move(breaks), a, b, opt);
}

NurbsSurface pipe_surface(const NurbsCurve& path, const std::function<double(double)>& radius,
                          const SweepOptions& opt) {
  const double a = path.domain_start(), b = path.domain_end();
  const auto r_at = [&](double v) {
    const double r = radius((v - a) / (b - a));
    if (!(r > 0.0)) throw GeometryError("pipe: the radius must be positive along the path");
    return r;
  };
  const Vec3 x0 = path.evaluate(a);
  const RotationMinimizingFrame rmf(path);
  const Mat3 F0 = rmf.frame(a);
  // unit circle at the path start, normal to the tangent; scaled to r(v) at every section
  const NurbsCurve profile = make_circle(x0, unit_tangent(path, a), 1.0);
  const Placement place = [&](double v) {
    return Transform::translation(path.evaluate(v)) * linear(relative(rmf.frame(v), F0)) *
           Transform::scaling({0, 0, 0}, r_at(v)) * Transform::translation(-x0);
  };
  return skin_adaptive(profile, place, break_points({&path}, a, b), opt);
}

}  // namespace nf
