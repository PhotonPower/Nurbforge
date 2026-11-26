#include "nurbforge/analysis/surface_analysis.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numeric>
#include <string>

#include "nurbforge/analysis/closest_point.hpp"
#include "nurbforge/analysis/properties.hpp"
#include "nurbforge/core/error.hpp"
#include "nurbforge/fitting/fit.hpp"

namespace nf {

namespace {

// samples_per_span + 1 parameters on every knot span of [a, b] (knots included once).
std::vector<double> grid(const KnotVector& U, double a, double b, std::size_t per_span,
                         const char* who) {
  if (per_span == 0) throw GeometryError(std::string(who) + ": samples_per_span must be > 0");
  std::vector<double> x{a};
  for (const double k : unique_knots(U).first)
    if (k > a + tol::kParametric && k < b - tol::kParametric) x.push_back(k);
  x.push_back(b);
  std::vector<double> out{a};
  for (std::size_t i = 0; i + 1 < x.size(); ++i)
    for (std::size_t q = 1; q <= per_span; ++q)
      out.push_back(q == per_span ? x[i + 1]
                                  : x[i] + (x[i + 1] - x[i]) * static_cast<double>(q) /
                                               static_cast<double>(per_span));
  return out;
}

struct RangeTracker {
  SurfaceMapRange r;
  bool any{false};
  void add(double value, double u, double v) {
    if (!any || value < r.min) {
      r.min = value;
      r.min_u = u;
      r.min_v = v;
    }
    if (!any || value > r.max) {
      r.max = value;
      r.max_u = u;
      r.max_v = v;
    }
    any = true;
  }
};

CurvatureMap curvature_map_of(const NurbsSurface& s, const TrimmedSurface* face,
                              std::size_t per_span) {
  const auto du = s.domain_u(), dv = s.domain_v();
  CurvatureMap m;
  m.u = grid(s.knots_u(), du[0], du[1], per_span, "curvature_map");
  m.v = grid(s.knots_v(), dv[0], dv[1], per_span, "curvature_map");
  RangeTracker K, H, k1, k2;
  double kmax = 0.0;
  for (const double u : m.u)
    for (const double v : m.v) {
      const SurfaceCurvature c = curvature(s, u, v);
      m.gaussian.push_back(c.gaussian);
      m.mean.push_back(c.mean);
      m.k_max.push_back(c.k1);
      m.k_min.push_back(c.k2);
      const bool in = face == nullptr || face->contains(u, v);
      m.inside.push_back(in);
      if (!in) continue;
      K.add(c.gaussian, u, v);
      H.add(c.mean, u, v);
      k1.add(c.k1, u, v);
      k2.add(c.k2, u, v);
      kmax = std::max({kmax, std::abs(c.k1), std::abs(c.k2)});
    }
  m.gaussian_range = K.r;
  m.mean_range = H.r;
  m.k_max_range = k1.r;
  m.k_min_range = k2.r;
  m.min_radius = kmax > 0.0 ? 1.0 / kmax : std::numeric_limits<double>::infinity();
  return m;
}

Vec3 unit_direction(const Vec3& d, const char* who) {
  const double n = norm(d);
  if (!(n > 0.0) || !std::isfinite(n))
    throw GeometryError(std::string(who) + ": the pull direction must be a non-zero vector");
  return d / n;
}

DraftAngleMap draft_map_of(const NurbsSurface& s, const TrimmedSurface* face, const Vec3& pull,
                           std::size_t per_span) {
  const Vec3 d = unit_direction(pull, "draft_angle_map");
  const auto du = s.domain_u(), dv = s.domain_v();
  DraftAngleMap m;
  m.u = grid(s.knots_u(), du[0], du[1], per_span, "draft_angle_map");
  m.v = grid(s.knots_v(), dv[0], dv[1], per_span, "draft_angle_map");
  RangeTracker r;
  for (const double u : m.u)
    for (const double v : m.v) {
      const double a = std::asin(std::clamp(dot(s.normal(u, v), d), -1.0, 1.0));
      m.angle.push_back(a);
      const bool in = face == nullptr || face->contains(u, v);
      m.inside.push_back(in);
      if (in) r.add(a, u, v);
    }
  m.range = r.r;
  return m;
}

// Curvature tensor k1 d1 d1^T + k2 d2 d2^T as a symmetric 3x3 matrix (6 entries), with the sign
// for the normal `n` (the normal curvature in a tangent direction t is t^T W t).
std::array<double, 6> curvature_tensor(const NurbsSurface& s, double u, double v, const Vec3& n) {
  const SurfaceCurvature c = curvature(s, u, v);
  const double sign = dot(c.normal, n) < 0.0 ? -1.0 : 1.0;
  const auto outer = [](const Vec3& a, double k) {
    return std::array<double, 6>{k * a.x * a.x, k * a.y * a.y, k * a.z * a.z,
                                 k * a.x * a.y, k * a.x * a.z, k * a.y * a.z};
  };
  const auto A = outer(c.dir1, sign * c.k1), B = outer(c.dir2, sign * c.k2);
  std::array<double, 6> W{};
  for (std::size_t i = 0; i < 6; ++i) W[i] = A[i] + B[i];
  return W;
}

double tensor_norm(const std::array<double, 6>& a, const std::array<double, 6>& b) {
  double s = 0.0;
  for (std::size_t i = 0; i < 6; ++i) {
    const double d = a[i] - b[i];
    s += (i < 3 ? 1.0 : 2.0) * d * d;  // Frobenius norm of the symmetric matrix
  }
  return std::sqrt(s);
}

// Accumulates the comparison of the two surfaces at a pair of corresponding points.
struct Joint {
  SurfaceContinuity r;
  double kmax{0.0};
  bool curvature_ok{true};

  void add(const NurbsSurface& a, double ua, double va, const NurbsSurface& b, double ub, double vb,
           double gap) {
    r.max_gap = std::max(r.max_gap, gap);
    try {
      const Vec3 na = a.normal(ua, va), nb = b.normal(ub, vb);
      r.max_angle = std::max(r.max_angle, std::asin(std::min(1.0, norm(cross(na, nb)))));
      const SurfaceCurvature ca = curvature(a, ua, va), cb = curvature(b, ub, vb);
      kmax = std::max({kmax, std::abs(ca.k1), std::abs(ca.k2), std::abs(cb.k1), std::abs(cb.k2)});
      // both tensors for the normal of a
      r.max_curvature_difference =
          std::max(r.max_curvature_difference,
                   tensor_norm(curvature_tensor(a, ua, va, na), curvature_tensor(b, ub, vb, na)));
    } catch (const GeometryError&) {
      curvature_ok = false;  // no normal or curvature at a degenerate point
    }
  }

  [[nodiscard]] SurfaceContinuity result(const SurfaceContinuityTolerance& tol) {
    r.geometric = -1;
    if (r.max_gap > tol.position) return r;
    r.geometric = 0;
    if (r.max_angle > tol.angle) return r;
    r.geometric = 1;
    if (curvature_ok && (r.max_curvature_difference <= tol.curvature * kmax || kmax <= 1e-12))
      r.geometric = 2;
    return r;
  }
};

}  // namespace

CurvatureMap curvature_map(const NurbsSurface& surface, std::size_t samples_per_span) {
  return curvature_map_of(surface, nullptr, samples_per_span);
}

CurvatureMap curvature_map(const TrimmedSurface& face, std::size_t samples_per_span) {
  return curvature_map_of(face.surface(), face.is_untrimmed() ? nullptr : &face, samples_per_span);
}

double draft_angle(const NurbsSurface& surface, double u, double v, const Vec3& pull_direction) {
  return std::asin(std::clamp(
      dot(surface.normal(u, v), unit_direction(pull_direction, "draft_angle")), -1.0, 1.0));
}

DraftAngleMap draft_angle_map(const NurbsSurface& surface, const Vec3& pull_direction,
                              std::size_t samples_per_span) {
  return draft_map_of(surface, nullptr, pull_direction, samples_per_span);
}

DraftAngleMap draft_angle_map(const TrimmedSurface& face, const Vec3& pull_direction,
                              std::size_t samples_per_span) {
  return draft_map_of(face.surface(), face.is_untrimmed() ? nullptr : &face, pull_direction,
                      samples_per_span);
}

SurfaceDirection surface_direction(const NurbsSurface& surface, double tolerance) {
  const auto du = surface.domain_u(), dv = surface.domain_v();
  const double um = 0.5 * (du[0] + du[1]), vm = 0.5 * (dv[0] + dv[1]);
  const auto D = surface.first_derivatives(um, vm);
  SurfaceDirection d;
  d.point = D[0];
  if (!(norm(D[1]) > 0.0) || !(norm(D[2]) > 0.0))
    throw GeometryError("surface_direction: the surface is degenerate at the middle of its domain");
  d.u_dir = normalized(D[1]);
  d.v_dir = normalized(D[2]);
  d.normal = surface.normal(um, vm);
  // sides compared at samples
  double close_u = 0.0, close_v = 0.0, ext[4] = {0.0, 0.0, 0.0, 0.0};
  const Vec3 c[4] = {surface.evaluate(du[0], dv[0]), surface.evaluate(du[1], dv[0]),
                     surface.evaluate(du[0], dv[0]), surface.evaluate(du[0], dv[1])};
  for (int i = 0; i <= 16; ++i) {
    const double u = du[0] + (du[1] - du[0]) * i / 16.0, v = dv[0] + (dv[1] - dv[0]) * i / 16.0;
    const Vec3 u0 = surface.evaluate(du[0], v), u1 = surface.evaluate(du[1], v);
    const Vec3 v0 = surface.evaluate(u, dv[0]), v1 = surface.evaluate(u, dv[1]);
    close_u = std::max(close_u, distance(u0, u1));
    close_v = std::max(close_v, distance(v0, v1));
    ext[0] = std::max(ext[0], distance(u0, c[0]));
    ext[1] = std::max(ext[1], distance(u1, c[1]));
    ext[2] = std::max(ext[2], distance(v0, c[2]));
    ext[3] = std::max(ext[3], distance(v1, c[3]));
  }
  for (int k = 0; k < 4; ++k) d.pole[k] = ext[k] <= tolerance;
  d.closed_u = close_u <= tolerance && !(d.pole[0] && d.pole[1]);
  d.closed_v = close_v <= tolerance && !(d.pole[2] && d.pole[3]);
  try {
    d.planar = fit_plane(surface.control_points()).max_deviation <= tolerance;
  } catch (const GeometryError&) {
    d.planar = true;  // collinear control points
  }
  return d;
}

SurfaceContinuity continuity(const NurbsSurface& a, SurfaceSide side, const NurbsSurface& b,
                             const SurfaceContinuityTolerance& tol, std::size_t samples) {
  if (samples < 2) throw GeometryError("continuity: at least two samples are needed");
  const auto du = a.domain_u(), dv = a.domain_v();
  Joint joint;
  for (std::size_t i = 0; i < samples; ++i) {
    // interior samples of the side (its ends are corners, where a third surface may meet)
    const double f = (static_cast<double>(i) + 0.5) / static_cast<double>(samples);
    double u = du[0] + f * (du[1] - du[0]), v = dv[0] + f * (dv[1] - dv[0]);
    if (side == SurfaceSide::UMin) u = du[0];
    if (side == SurfaceSide::UMax) u = du[1];
    if (side == SurfaceSide::VMin) v = dv[0];
    if (side == SurfaceSide::VMax) v = dv[1];
    const SurfacePoint q = closest_point(b, a.evaluate(u, v));
    joint.add(a, u, v, b, q.u, q.v, q.distance);
  }
  return joint.result(tol);
}

std::vector<SurfaceContinuity> edge_continuity(const Shell& shell,
                                               const SurfaceContinuityTolerance& tol,
                                               std::size_t samples) {
  if (samples == 0) throw GeometryError("edge_continuity: at least one sample is needed");
  std::vector<SurfaceContinuity> out;
  out.reserve(shell.edges().size());
  for (const ShellEdge& e : shell.edges()) {
    if (e.uses.size() != 2) {
      out.push_back({});
      continue;
    }
    const double a = e.curve.domain_start(), b = e.curve.domain_end();
    Joint joint;
    for (std::size_t i = 0; i < samples; ++i) {
      const double fr = (static_cast<double>(i) + 0.5) / static_cast<double>(samples);
      const Vec3 P = e.curve.evaluate(a + fr * (b - a));
      SurfacePoint sp[2];
      const NurbsSurface* S[2] = {nullptr, nullptr};
      for (std::size_t k = 0; k < 2; ++k) {
        const Coedge& use = e.uses[k];
        const TrimmedSurface& f = shell.faces()[use.face];
        const TrimLoop& loop = use.loop == 0 ? f.outer() : f.holes()[use.loop - 1];
        const NurbsCurve& pc = loop.curves[use.index];
        const double tau = use.reversed ? 1.0 - fr : fr;
        const Vec3 uv =
            pc.evaluate(pc.domain_start() + tau * (pc.domain_end() - pc.domain_start()));
        S[k] = &f.surface();
        sp[k] = closest_point_local(f.surface(), P, uv.x, uv.y);
      }
      joint.add(*S[0], sp[0].u, sp[0].v, *S[1], sp[1].u, sp[1].v,
                distance(sp[0].point, sp[1].point));
    }
    out.push_back(joint.result(tol));
  }
  return out;
}

std::vector<std::vector<std::size_t>> find_duplicate_surfaces(
    const std::vector<NurbsSurface>& surfaces, double tolerance) {
  if (!(tolerance > 0.0))
    throw GeometryError("find_duplicate_surfaces: tolerance must be positive");
  const std::size_t n = surfaces.size();
  struct Info {
    Vec3 lo, hi;
    std::vector<Vec3> samples;
  };
  std::vector<Info> info;
  constexpr int kN = 6;
  for (const auto& s : surfaces) {
    Info f{{1e300, 1e300, 1e300}, {-1e300, -1e300, -1e300}, {}};
    for (const auto& p : s.control_points()) {
      f.lo = {std::min(f.lo.x, p.x), std::min(f.lo.y, p.y), std::min(f.lo.z, p.z)};
      f.hi = {std::max(f.hi.x, p.x), std::max(f.hi.y, p.y), std::max(f.hi.z, p.z)};
    }
    const auto du = s.domain_u(), dv = s.domain_v();
    for (int i = 0; i <= kN; ++i)
      for (int j = 0; j <= kN; ++j)
        f.samples.push_back(
            s.evaluate(du[0] + (du[1] - du[0]) * i / kN, dv[0] + (dv[1] - dv[0]) * j / kN));
    info.push_back(std::move(f));
  }
  const auto on = [&](std::size_t i, std::size_t j) {
    for (const auto& p : info[i].samples)
      if (closest_point(surfaces[j], p).distance > tolerance) return false;
    return true;
  };
  const auto same = [&](std::size_t i, std::size_t j) {
    const Info& a = info[i];
    const Info& b = info[j];
    if (a.lo.x > b.hi.x + tolerance || b.lo.x > a.hi.x + tolerance || a.lo.y > b.hi.y + tolerance ||
        b.lo.y > a.hi.y + tolerance || a.lo.z > b.hi.z + tolerance || b.lo.z > a.hi.z + tolerance)
      return false;
    return on(i, j) && on(j, i);
  };
  std::vector<std::size_t> parent(n);
  std::iota(parent.begin(), parent.end(), std::size_t{0});
  const auto find = [&](std::size_t x) {
    while (parent[x] != x) x = parent[x] = parent[parent[x]];
    return x;
  };
  for (std::size_t i = 0; i < n; ++i)
    for (std::size_t j = i + 1; j < n; ++j) {
      const std::size_t ri = find(i), rj = find(j);
      if (ri != rj && same(i, j)) parent[std::max(ri, rj)] = std::min(ri, rj);
    }
  std::vector<std::vector<std::size_t>> groups(n);
  for (std::size_t i = 0; i < n; ++i) groups[find(i)].push_back(i);
  std::vector<std::vector<std::size_t>> out;
  for (auto& g : groups)
    if (g.size() > 1) out.push_back(std::move(g));
  return out;
}

}  // namespace nf
