#include "nurbforge/modeling/blend.hpp"

#include <cmath>
#include <vector>

#include "nurbforge/core/error.hpp"
#include "nurbforge/nurbs/knot_vector.hpp"

namespace nf {

namespace {

struct EndData {
  Vec3 p;  // end point
  Vec3 t;  // unit tangent, pointing along the blend direction
  Vec3 k;  // curvature vector (independent of the direction of travel)
};

EndData end_data(const NurbsCurve& c, CurveEnd e, bool leaving) {
  const double u = e == CurveEnd::End ? c.domain_end() : c.domain_start();
  const auto d = c.derivatives(u, 2);
  const double s = norm(d[1]);
  if (!(s > 0.0)) throw GeometryError("blend_curves: vanishing tangent at a curve end");
  EndData r;
  r.p = d[0];
  // a: leave the curve (forward at its end, backward at its start);
  // b: enter the curve towards its interior (forward at its start, backward at its end)
  const bool forward = leaving ? e == CurveEnd::End : e == CurveEnd::Start;
  r.t = (forward ? 1.0 : -1.0) / s * d[1];
  const Vec3 c2 = cross(d[1], d[2]);
  r.k = (1.0 / (s * s * s * s)) * cross(c2, d[1]);  // (C' x C'') x C' / |C'|^4
  return r;
}

}  // namespace

NurbsCurve blend_curves(const NurbsCurve& a, CurveEnd end0, const NurbsCurve& b, CurveEnd end1,
                        int continuity, double tension) {
  if (continuity != 1 && continuity != 2)
    throw GeometryError("blend_curves: continuity must be 1 (G1) or 2 (G2)");
  if (!(tension > 0.0)) throw GeometryError("blend_curves: tension must be positive");
  const EndData s = end_data(a, end0, true), e = end_data(b, end1, false);
  const double chord = distance(s.p, e.p);
  if (!(chord > 0.0)) throw GeometryError("blend_curves: the curve ends coincide");
  std::vector<Vec3> P;
  if (continuity == 1) {
    // cubic: C'(0) = 3 (P1 - P0) = m t0, m = chord
    const double m = tension * chord;
    P = {s.p, s.p + (m / 3.0) * s.t, e.p - (m / 3.0) * e.t, e.p};
  } else {
    // quintic: C'(0) = 5 (P1 - P0) = m t0, C''(0) = 20 (P2 - 2 P1 + P0) = m^2 k0 (G2 with
    // zero tangential acceleration); symmetric at the end
    const double m = tension * chord;
    const Vec3 p1 = s.p + (m / 5.0) * s.t, q1 = e.p - (m / 5.0) * e.t;
    const Vec3 p2 = 2.0 * p1 - s.p + (m * m / 20.0) * s.k;
    const Vec3 q2 = 2.0 * q1 - e.p + (m * m / 20.0) * e.k;
    P = {s.p, p1, p2, q2, q1, e.p};
  }
  const std::size_t p = P.size() - 1;
  return {p, clamped_uniform_knots(p, P.size()), P, {}};
}

}  // namespace nf
