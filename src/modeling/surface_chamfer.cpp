#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

#include "fit_family.hpp"
#include "geodesic.hpp"
#include "nurbforge/analysis/closest_point.hpp"
#include "nurbforge/analysis/surface_intersect.hpp"
#include "nurbforge/core/error.hpp"
#include "nurbforge/modeling/construct.hpp"
#include "nurbforge/modeling/fillet.hpp"

namespace nf {

SurfaceChamfer chamfer_surfaces(const NurbsSurface& s1, const NurbsSurface& s2, double distance1,
                                double distance2, const SurfaceChamferOptions& opt) {
  if (!(distance1 > 0.0) || !(distance2 > 0.0))
    throw GeometryError("chamfer_surfaces: the distances must be positive");
  if ((opt.side1 != 1 && opt.side1 != -1) || (opt.side2 != 1 && opt.side2 != -1))
    throw GeometryError("chamfer_surfaces: side1 and side2 must be +1 or -1");
  if (!(opt.tolerance > 0.0)) throw GeometryError("chamfer_surfaces: tolerance must be positive");

  const SurfaceIntersection cut = intersect(s1, s2);
  if (cut.curves.empty()) throw GeometryError("chamfer_surfaces: the surfaces do not intersect");
  const IntersectionCurve* edge = &cut.curves.front();
  const auto length = [](const IntersectionCurve& c) {
    double l = 0.0;
    for (std::size_t i = 1; i < c.points.size(); ++i)
      l += distance(c.points[i - 1].point, c.points[i].point);
    return l;
  };
  const auto gap = [&](const IntersectionCurve& c) {
    double d = std::numeric_limits<double>::infinity();
    for (const auto& p : c.points) d = std::min(d, distance(p.point, *opt.near));
    return d;
  };
  for (const auto& c : cut.curves)
    if (opt.near ? gap(c) < gap(*edge) : length(c) > length(*edge)) edge = &c;

  const double a = edge->curve.domain_start(), b = edge->curve.domain_end();
  int steps1 = 8, steps2 = 8;
  // the contact points and their parameters for the edge parameter t
  const auto exact = [&](double t) {
    const Vec3 q1 = edge->pcurve1.evaluate(t), q2 = edge->pcurve2.evaluate(t);
    // the exact edge point (the fitted intersection curve is within its own tolerance only)
    SurfacePoint p1 = closest_point_local(s1, edge->curve.evaluate(t), q1.x, q1.y);
    SurfacePoint p2 = closest_point_local(s2, p1.point, q2.x, q2.y);
    for (int it = 0; it < 100 && distance(p1.point, p2.point) > 1e-14 * (1.0 + norm(p1.point));
         ++it) {
      p1 = closest_point_local(s1, p2.point, p1.u, p1.v);
      p2 = closest_point_local(s2, p1.point, p2.u, p2.v);
    }
    const Vec3 n1 = s1.normal(p1.u, p1.v), n2 = s2.normal(p2.u, p2.v);
    const Vec3 T = cross(n1, n2);
    if (!(norm(T) > 1e-9))
      throw GeometryError("chamfer_surfaces: the surfaces are tangent at t = " + std::to_string(t));
    // directions within the surfaces, normal to the edge, into the wedge of the chosen sides
    Vec3 d1 = cross(n1, T), d2 = cross(n2, T);
    if (dot(d1, n2) * opt.side2 < 0.0) d1 = -d1;
    if (dot(d2, n1) * opt.side1 < 0.0) d2 = -d2;
    const auto w1 =
        detail::geodesic_point(s1, p1.u, p1.v, d1, distance1, 0.5 * opt.tolerance, steps1);
    const auto w2 =
        detail::geodesic_point(s2, p2.u, p2.v, d2, distance2, 0.5 * opt.tolerance, steps2);
    return std::vector<Vec3>{s1.evaluate(w1[0], w1[1]), s2.evaluate(w2[0], w2[1]),
                             Vec3{w1[0], w1[1], 0.0}, Vec3{w2[0], w2[1], 0.0}};
  };
  const auto error = [&](const std::vector<NurbsCurve>& c, double t, const std::vector<Vec3>& e) {
    const Vec3 r1 = c[2].evaluate(t), r2 = c[3].evaluate(t);
    const auto d1 = s1.domain_u(), e1 = s1.domain_v(), d2 = s2.domain_u(), e2 = s2.domain_v();
    return std::max(
        {distance(c[0].evaluate(t), e[0]), distance(c[1].evaluate(t), e[1]),
         distance(s1.evaluate(std::clamp(r1.x, d1[0], d1[1]), std::clamp(r1.y, e1[0], e1[1])),
                  e[0]),
         distance(s2.evaluate(std::clamp(r2.x, d2[0], d2[1]), std::clamp(r2.y, e2[0], e2[1])),
                  e[1])});
  };
  // start with four samples per knot span of the edge (tiny spans of the fitted curve merged)
  std::vector<double> T{a};
  const auto knots = unique_knots(edge->curve.knots()).first;
  for (std::size_t i = 0; i + 1 < knots.size(); ++i) {
    if (knots[i + 1] <= a || knots[i] >= b || knots[i + 1] - knots[i] < 1e-3 * (b - a)) continue;
    for (int j = 1; j <= 4; ++j) {
      const double t = knots[i] + (knots[i + 1] - knots[i]) * j / 4.0;
      if (t > T.back() + 1e-9 * (b - a) && t < b - 1e-9 * (b - a)) T.push_back(t);
    }
  }
  T.push_back(b);
  auto c = detail::fit_family(T, exact, error, 0.5 * opt.tolerance, "chamfer_surfaces");
  NurbsSurface surface = ruled(c[0], c[1]).transposed();  // u across, v along the edge
  const auto dv = surface.domain_v();
  for (auto& x : c) x = reparametrize(x, dv[0], dv[1]);
  return {std::move(surface), reparametrize(edge->curve, dv[0], dv[1]),
          std::move(c[0]),    std::move(c[1]),
          std::move(c[2]),    std::move(c[3])};
}

}  // namespace nf
