#include "nurbforge/analysis/curve_analysis.hpp"

#include <algorithm>
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

double size_of(const NurbsCurve& c) {
  double size = 0.0;
  for (const auto& p : c.control_points())
    size = std::max(size, distance(p, c.control_points().front()));
  return size;
}

// Point, unit tangent and curvature vector from the derivatives {C, C', C''}.
CurveCurvature from_derivatives(const std::vector<Vec3>& d, double u) {
  const double speed = norm(d[1]);
  if (!(speed > 0.0))
    throw GeometryError("curve_curvature: the derivative vanishes at u = " + std::to_string(u));
  CurveCurvature r;
  r.point = d[0];
  r.tangent = d[1] / speed;
  r.vector = (d[2] - dot(d[2], r.tangent) * r.tangent) / (speed * speed);
  r.curvature = norm(r.vector);
  r.center = r.point;
  if (r.curvature > 1e-14 * std::max(1.0, 1.0 / std::max(speed, 1e-300))) {
    r.normal = r.vector / r.curvature;
    r.radius = 1.0 / r.curvature;
    r.center = r.point + r.radius * r.normal;
  } else {
    r.curvature = 0.0;
    r.vector = {0, 0, 0};
    r.radius = std::numeric_limits<double>::infinity();
  }
  return r;
}

}  // namespace

CurveCurvature curve_curvature(const NurbsCurve& c, double u) {
  return from_derivatives(c.derivatives(u, 2), u);
}

CurvatureGraph curvature_graph(const NurbsCurve& c, double scale, std::size_t samples_per_span) {
  if (samples_per_span == 0) throw GeometryError("curvature_graph: samples_per_span must be > 0");
  CurvatureGraph g;
  for (const NurbsCurve& seg : c.decompose_bezier()) {  // one-sided values at the knots
    const double a = seg.domain_start(), b = seg.domain_end();
    for (std::size_t k = 0; k <= samples_per_span; ++k) {
      const double u = k == samples_per_span ? b
                                             : a + (b - a) * static_cast<double>(k) /
                                                       static_cast<double>(samples_per_span);
      const CurveCurvature cc = from_derivatives(seg.derivatives(u, 2), u);
      g.parameters.push_back(u);
      g.points.push_back(cc.point);
      g.tips.push_back(cc.point - scale * cc.vector);
      g.curvatures.push_back(cc.curvature);
      if (cc.curvature > g.max_curvature) {
        g.max_curvature = cc.curvature;
        g.max_parameter = u;
      }
    }
  }
  return g;
}

bool is_closed(const NurbsCurve& c, double tolerance) {
  return distance(c.evaluate(c.domain_start()), c.evaluate(c.domain_end())) <= tolerance;
}

bool is_planar(const NurbsCurve& c, double tolerance) {
  if (c.control_points().size() <= 3) {
    // up to three points always lie in a plane
    return true;
  }
  try {
    return fit_plane(c.control_points()).max_deviation <= tolerance;
  } catch (const GeometryError&) {
    return true;  // collinear or coincident control points
  }
}

CurveDirection curve_direction(const NurbsCurve& c, const Vec3& plane_normal, double tolerance) {
  const double a = c.domain_start(), b = c.domain_end();
  CurveDirection d;
  const auto da = c.derivatives(a, 1), db = c.derivatives(b, 1);
  d.start = da[0];
  d.end = db[0];
  if (!(norm(da[1]) > 0.0) || !(norm(db[1]) > 0.0))
    throw GeometryError("curve_direction: the derivative vanishes at an end of the curve");
  d.start_tangent = normalized(da[1]);
  d.end_tangent = normalized(db[1]);
  d.closed = is_closed(c, tolerance);
  d.planar = is_planar(c, tolerance);
  if (norm(plane_normal) > 0.0) {
    d.plane_normal = normalized(plane_normal);
  } else if (d.planar) {
    try {
      d.plane_normal = normalized(fit_plane(c.control_points()).plane.normal);
    } catch (const GeometryError&) {
      // coincident control points: no plane
    }
  }
  if (d.closed && norm(d.plane_normal) > 0.0) {
    // enclosed area of the projection: 1/2 int (C - O) x C' . n
    const Vec3 origin = d.start, n = d.plane_normal;
    std::vector<double> breaks;
    for (const double k : unique_knots(c.knots()).first)
      if (k > a && k < b) breaks.push_back(k);
    d.signed_area = 0.5 * integrate(
                              [&](double u) {
                                const auto e = c.derivatives(u, 1);
                                return dot(cross(e[0] - origin, e[1]), n);
                              },
                              a, b, 1e-10, breaks);
    const double size = size_of(c);
    if (std::abs(d.signed_area) > 1e-12 * size * size) d.orientation = d.signed_area > 0.0 ? 1 : -1;
  }
  return d;
}

namespace {

// {C, C', C''} at an end, in the direction of travel through a joint: `forward` if the curve
// is traversed as parametrised.
std::vector<Vec3> end_derivatives(const NurbsCurve& c, CurveEnd end, bool forward) {
  auto d = c.derivatives(end == CurveEnd::Start ? c.domain_start() : c.domain_end(), 2);
  if (!forward) d[1] = -d[1];  // reversed parameter: odd derivatives change sign
  return d;
}

CurveContinuity compare(const std::vector<Vec3>& da, const std::vector<Vec3>& db,
                        const ContinuityTolerance& tol) {
  CurveContinuity r;
  r.gap = distance(da[0], db[0]);
  const CurveCurvature ca = from_derivatives(da, 0.0), cb = from_derivatives(db, 0.0);
  r.tangent_angle = std::atan2(norm(cross(ca.tangent, cb.tangent)), dot(ca.tangent, cb.tangent));
  r.curvature_a = ca.curvature;
  r.curvature_b = cb.curvature;
  r.curvature_difference = norm(ca.vector - cb.vector);
  if (r.gap > tol.position) return r;
  r.geometric = r.parametric = 0;
  if (r.tangent_angle <= tol.angle) {
    r.geometric = 1;
    const double kmax = std::max(ca.curvature, cb.curvature);
    if (r.curvature_difference <= tol.curvature * kmax || kmax <= 1e-12) r.geometric = 2;
  }
  // parametric: the derivative vectors agree, relative to their size
  const double s1 = std::max(norm(da[1]), norm(db[1]));
  if (norm(da[1] - db[1]) <= std::max(tol.angle, 1e-9) * s1) {
    r.parametric = 1;
    const double s2 = std::max({norm(da[2]), norm(db[2]), s1});
    if (norm(da[2] - db[2]) <= std::max(tol.curvature, 1e-9) * s2) r.parametric = 2;
  }
  return r;
}

}  // namespace

CurveContinuity continuity(const NurbsCurve& a, CurveEnd end_a, const NurbsCurve& b, CurveEnd end_b,
                           const ContinuityTolerance& tol) {
  // travel from a into b: a forwards if it ends at the joint, b forwards if it starts there
  return compare(end_derivatives(a, end_a, end_a == CurveEnd::End),
                 end_derivatives(b, end_b, end_b == CurveEnd::Start), tol);
}

CurveContinuity continuity(const NurbsCurve& a, const NurbsCurve& b,
                           const ContinuityTolerance& tol) {
  const auto at = [](const NurbsCurve& c, CurveEnd e) {
    return c.evaluate(e == CurveEnd::Start ? c.domain_start() : c.domain_end());
  };
  CurveEnd ea = CurveEnd::End, eb = CurveEnd::Start;
  double best = std::numeric_limits<double>::infinity();
  // the usual chain (end of a to start of b) wins ties
  for (const CurveEnd x : {CurveEnd::End, CurveEnd::Start})
    for (const CurveEnd y : {CurveEnd::Start, CurveEnd::End}) {
      const double gap = distance(at(a, x), at(b, y));
      if (gap < best - 1e-14 * (1.0 + best)) {
        best = gap;
        ea = x;
        eb = y;
      }
    }
  return continuity(a, ea, b, eb, tol);
}

std::vector<CurveBreak> continuity_breaks(const NurbsCurve& c, const ContinuityTolerance& tol) {
  std::vector<CurveBreak> out;
  const auto segs = c.decompose_bezier();
  const auto smooth = [](const CurveContinuity& k) { return k.geometric == 2; };
  for (std::size_t i = 0; i + 1 < segs.size(); ++i) {
    const CurveContinuity k = continuity(segs[i], CurveEnd::End, segs[i + 1], CurveEnd::Start, tol);
    if (!smooth(k)) out.push_back({segs[i].domain_end(), k});
  }
  if (is_closed(c, tol.position) && !segs.empty()) {
    const CurveContinuity k =
        continuity(segs.back(), CurveEnd::End, segs.front(), CurveEnd::Start, tol);
    if (!smooth(k)) out.push_back({c.domain_end(), k});
  }
  return out;
}

std::vector<std::vector<std::size_t>> find_duplicate_curves(const std::vector<NurbsCurve>& curves,
                                                            double tolerance) {
  if (!(tolerance > 0.0)) throw GeometryError("find_duplicate_curves: tolerance must be positive");
  const std::size_t n = curves.size();
  struct Info {
    Vec3 lo, hi;
    double length;
    std::vector<Vec3> samples;
  };
  std::vector<Info> info;
  info.reserve(n);
  constexpr int kSamples = 24;
  for (const auto& c : curves) {
    Info f{{1e300, 1e300, 1e300}, {-1e300, -1e300, -1e300}, arc_length(c), {}};
    for (const auto& p : c.control_points()) {
      f.lo = {std::min(f.lo.x, p.x), std::min(f.lo.y, p.y), std::min(f.lo.z, p.z)};
      f.hi = {std::max(f.hi.x, p.x), std::max(f.hi.y, p.y), std::max(f.hi.z, p.z)};
    }
    const double a = c.domain_start(), b = c.domain_end();
    for (int k = 0; k <= kSamples; ++k) f.samples.push_back(c.evaluate(a + (b - a) * k / kSamples));
    info.push_back(std::move(f));
  }
  const auto on = [&](std::size_t i, std::size_t j) {  // the samples of i lie on j
    for (const auto& p : info[i].samples)
      if (closest_point(curves[j], p).distance > tolerance) return false;
    return true;
  };
  const auto same = [&](std::size_t i, std::size_t j) {
    const Info& a = info[i];
    const Info& b = info[j];
    // equal lengths (each sample may be off by the tolerance) and overlapping boxes first
    if (std::abs(a.length - b.length) > 2.0 * tolerance * (kSamples + 1)) return false;
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
