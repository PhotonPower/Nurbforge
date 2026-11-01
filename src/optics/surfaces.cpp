#include "nurbforge/optics/surfaces.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numbers>
#include <optional>
#include <string>

#include "nurbforge/core/error.hpp"
#include "nurbforge/fitting/fit.hpp"
#include "nurbforge/geometry/primitives.hpp"
#include "nurbforge/modeling/construct.hpp"
#include "nurbforge/nurbs/join.hpp"

namespace nf {

// ------------------------------------------------------------------ sag functions

double Asphere::sag(double r) const {
  const double r2 = r * r, c = curvature;
  const double arg = 1.0 - (1.0 + conic) * c * c * r2;
  if (arg < 0.0) throw GeometryError("asphere: r = " + std::to_string(r) + " is beyond the conic's domain");
  double z = c * r2 / (1.0 + std::sqrt(arg));
  double rp = r2 * r2;
  for (double a : even) {
    z += a * rp;
    rp *= r2;
  }
  const double ar = std::abs(r);
  rp = ar;
  for (double b : odd) {
    z += b * rp;
    rp *= ar;
  }
  return z;
}

double Asphere::max_radius() const {
  if (curvature == 0.0 || 1.0 + conic <= 0.0) return std::numeric_limits<double>::infinity();
  return 1.0 / (std::abs(curvature) * std::sqrt(1.0 + conic));
}

void zernike_nm(std::size_t j, ZernikeIndexing indexing, int& n, int& m) {
  if (j < 1) throw GeometryError("zernike: index must be >= 1");
  if (indexing == ZernikeIndexing::Noll) {
    std::size_t nn = 0, j1 = j - 1;
    while (j1 > nn) {
      ++nn;
      j1 -= nn;
    }
    const std::size_t m_abs = (nn % 2) + 2 * ((j1 + ((nn + 1) % 2)) / 2);
    n = static_cast<int>(nn);
    m = m_abs == 0 ? 0 : (j % 2 == 0 ? static_cast<int>(m_abs) : -static_cast<int>(m_abs));
    return;
  }
  // Fringe: group d = (n + |m|) / 2 holds j = d^2 + 1 .. (d+1)^2, ordered m = d, -d, d-1, ..., 1, -1, 0
  std::size_t d = 0;
  while ((d + 1) * (d + 1) < j) ++d;
  const std::size_t pos = j - d * d - 1;  // 0 .. 2d
  const std::size_t m_abs = d - pos / 2;
  m = m_abs == 0 ? 0 : (pos % 2 == 0 ? static_cast<int>(m_abs) : -static_cast<int>(m_abs));
  n = static_cast<int>(2 * d - m_abs);
}

double zernike(std::size_t j, ZernikeIndexing indexing, double rho, double theta) {
  int n = 0, m = 0;
  zernike_nm(j, indexing, n, m);
  const int ma = std::abs(m);
  double R = 0.0;
  for (int s = 0; s <= (n - ma) / 2; ++s) {
    const double coef = std::tgamma(n - s + 1.0) /
                        (std::tgamma(s + 1.0) * std::tgamma((n + ma) / 2 - s + 1.0) *
                         std::tgamma((n - ma) / 2 - s + 1.0));
    R += (s % 2 == 0 ? coef : -coef) * std::pow(rho, n - 2 * s);
  }
  double ang = 1.0;
  if (m > 0) ang = std::cos(m * theta);
  if (m < 0) ang = std::sin(ma * theta);
  double norm = 1.0;
  if (indexing == ZernikeIndexing::Noll) norm = std::sqrt((m == 0 ? 1.0 : 2.0) * (n + 1));
  return norm * R * ang;
}

double ZernikeSag::sag(double x, double y) const {
  const double r = std::hypot(x, y);
  double z = base.sag(r);
  const double rho = r / norm_radius, theta = std::atan2(y, x);
  for (std::size_t j = 0; j < coeffs.size(); ++j)
    if (coeffs[j] != 0.0) z += coeffs[j] * zernike(j + 1, indexing, rho, theta);
  return z;
}

double XYPolynomialSag::sag(double x, double y) const {
  double z = base.sag(std::hypot(x, y));
  const double xn = x / norm_radius, yn = y / norm_radius;
  for (const auto& t : terms) z += t.c * std::pow(xn, t.i) * std::pow(yn, t.j);
  return z;
}

namespace {

constexpr double pi = std::numbers::pi;

// ------------------------------------------------------------------ exact conic profile

NurbsCurve conic_profile(double c, double k, double r_max) {
  if (c == 0.0) return make_line({0, 0, 0}, {r_max, 0, 0});
  const auto z = [&](double r) { return c * r * r / (1.0 + std::sqrt(std::max(0.0, 1.0 - (1.0 + k) * c * c * r * r))); };
  const auto tangent = [&](double r) {  // (dr, dz), finite also where the tangent is vertical
    return std::array<double, 2>{std::sqrt(std::max(0.0, 1.0 - (1.0 + k) * c * c * r * r)), c * r};
  };
  const auto t_end = tangent(r_max);
  const double phi_max = std::atan2(std::abs(t_end[1]), t_end[0]);
  const auto n = static_cast<std::size_t>(std::max(1.0, std::ceil(phi_max / (pi / 4) - 1e-12)));
  std::vector<double> rb(n + 1);
  for (std::size_t i = 0; i <= n; ++i) {
    const double phi = phi_max * static_cast<double>(i) / static_cast<double>(n);
    const double s = std::sin(phi), co = std::cos(phi);
    rb[i] = s / (std::abs(c) * std::sqrt(co * co + (1.0 + k) * s * s));
  }
  rb[0] = 0.0;
  rb[n] = r_max;
  const auto G = [&](double r, double zz) { return c * r * r - 2.0 * zz + (1.0 + k) * c * zz * zz; };

  std::optional<NurbsCurve> profile;
  for (std::size_t i = 0; i < n; ++i) {
    const double ra = rb[i], rbb = rb[i + 1];
    const double x0 = ra, z0 = z(ra), x2 = rbb, z2 = z(rbb);
    const auto T0 = tangent(ra), T2 = tangent(rbb);
    // P1: intersection of the end tangents
    const double det = T0[0] * (-T2[1]) - T0[1] * (-T2[0]);
    const double a = ((x2 - x0) * (-T2[1]) - (z2 - z0) * (-T2[0])) / det;
    const double x1 = x0 + a * T0[0], z1 = z0 + a * T0[1];
    // weight from the shoulder point: S = M + s (P1 - M) on the conic, w = s / (1 - s)
    const double mx = 0.5 * (x0 + x2), mz = 0.5 * (z0 + z2);
    const auto g = [&](double s) { return G(mx + s * (x1 - mx), mz + s * (z1 - mz)); };
    const double g0 = g(0.0), g1 = g(0.5), g2 = g(1.0);
    const double A = 2.0 * (g0 - 2.0 * g1 + g2), B = -3.0 * g0 + 4.0 * g1 - g2, C = g0;
    double s = -1.0;
    if (std::abs(A) < 1e-300) {
      s = -C / B;
    } else {
      const double disc = std::sqrt(std::max(0.0, B * B - 4.0 * A * C));
      for (double root : {(-B + disc) / (2.0 * A), (-B - disc) / (2.0 * A)})
        if (root > 0.0 && root < 1.0) s = root;
    }
    if (!(s > 0.0 && s < 1.0)) throw GeometryError("asphere: cannot construct the conic profile");
    const double w = s / (1.0 - s);
    const double u0 = static_cast<double>(i), u1 = u0 + 1.0;
    const NurbsCurve seg{2, {u0, u0, u0, u1, u1, u1}, {{x0, 0, z0}, {x1, 0, z1}, {x2, 0, z2}}, {1.0, w, 1.0}};
    profile = profile ? join(*profile, seg, 1e-12) : seg;
  }
  return *profile;
}

// ------------------------------------------------------------------ graph fits

KnotVector knots_from_breaks(const std::vector<double>& br, std::size_t p) {
  KnotVector U(p + 1, br.front());
  for (std::size_t i = 1; i + 1 < br.size(); ++i) U.push_back(br[i]);
  U.insert(U.end(), p + 1, br.back());
  return U;
}

// Fit parameters: p + 2 points per span (shared ends).
std::vector<double> fit_params(const std::vector<double>& br, std::size_t p) {
  std::vector<double> t{br.front()};
  const std::size_t per = p + 2;
  for (std::size_t i = 0; i + 1 < br.size(); ++i)
    for (std::size_t q = 1; q <= per; ++q)
      t.push_back(q == per ? br[i + 1] : br[i] + (br[i + 1] - br[i]) * static_cast<double>(q) / static_cast<double>(per));
  return t;
}

void refine(std::vector<double>& br, const std::vector<bool>& bad) {
  std::vector<double> out{br.front()};
  for (std::size_t i = 0; i + 1 < br.size(); ++i) {
    if (bad[i]) out.push_back(0.5 * (br[i] + br[i + 1]));
    out.push_back(br[i + 1]);
  }
  br = std::move(out);
}

double checked(double v) {
  if (!std::isfinite(v)) throw GeometryError("sag function is not finite in the fitted region");
  return v;
}

// Cubic graph z(r) on [0, r_max]; symmetric fit on [-r_max, r_max] (zero slope at the vertex)
// unless `symmetric` is false.
NurbsCurve fit_profile(const std::function<double(double)>& f, double r_max, double tol, bool symmetric) {
  const std::size_t p = 3;
  std::vector<double> half{0.0, 0.5 * r_max, r_max};
  for (int iter = 0; iter < 40; ++iter) {
    std::vector<double> br;
    if (symmetric)
      for (std::size_t i = half.size(); i-- > 1;) br.push_back(-half[i]);
    br.insert(br.end(), half.begin(), half.end());
    const auto U = knots_from_breaks(br, p);
    const auto t = fit_params(br, p);
    std::vector<Vec3> pts;
    for (double x : t) pts.push_back({x, 0.0, checked(f(std::abs(x)))});
    const NurbsCurve c = approximate_curve(pts, p, U, t);
    std::vector<bool> bad(half.size() - 1, false);
    bool any = false;
    for (std::size_t i = 0; i + 1 < half.size(); ++i)
      for (int q = 0; q <= 8; ++q) {
        const double x = half[i] + (half[i + 1] - half[i]) * q / 8.0;
        const Vec3 s = c.evaluate(x);
        if (std::abs(s.z - f(x)) > tol || std::abs(s.x - x) > tol) {
          bad[i] = true;
          any = true;
          break;
        }
      }
    if (!any) return symmetric ? c.sub_curve(0.0, r_max) : c;
    if (half.size() > 4096) break;
    refine(half, bad);
  }
  throw GeometryError("asphere: profile tolerance not reached (too many spans)");
}

}  // namespace

NurbsCurve asphere_profile(const Asphere& a, double r_max, double tolerance) {
  if (!(r_max > 0.0)) throw GeometryError("asphere: r_max must be positive");
  if (r_max > a.max_radius() * (1.0 + 1e-12))
    throw GeometryError("asphere: r_max exceeds the conic's domain (" + std::to_string(a.max_radius()) + ")");
  r_max = std::min(r_max, a.max_radius());
  if (a.is_conic()) return conic_profile(a.curvature, a.conic, r_max);
  // symmetric fit (zero slope at the vertex) unless a linear term B1*r makes the vertex a cone tip
  const bool symmetric = a.odd.empty() || a.odd.front() == 0.0;
  return fit_profile([&](double r) { return a.sag(r); }, r_max, tolerance, symmetric);
}

NurbsSurface asphere_surface(const Asphere& a, double r_max, double tolerance) {
  return revolve(asphere_profile(a, r_max, tolerance), {0, 0, 0}, {0, 0, 1}, 2 * pi);
}

NurbsSurface sag_surface(const std::function<double(double, double)>& f, double x0, double x1,
                         double y0, double y1, std::size_t p, double tol, std::size_t max_spans) {
  if (!(x1 > x0) || !(y1 > y0)) throw GeometryError("sag_surface: empty rectangle");
  if (p < 1) throw GeometryError("sag_surface: degree must be >= 1");
  std::vector<double> bx{x0, 0.5 * (x0 + x1), x1}, by{y0, 0.5 * (y0 + y1), y1};
  while (true) {
    const auto U = knots_from_breaks(bx, p), V = knots_from_breaks(by, p);
    const auto xs = fit_params(bx, p), ys = fit_params(by, p);
    // rows (fixed y) along x, then columns of the intermediate control points along y
    std::vector<std::vector<Vec3>> R;  // R[j][i]
    for (double y : ys) {
      std::vector<Vec3> row;
      for (double x : xs) row.push_back({x, y, checked(f(x, y))});
      R.push_back(approximate_curve(row, p, U, xs).control_points());
    }
    const std::size_t nx = R.front().size();
    std::vector<std::vector<Vec3>> C(nx);
    for (std::size_t i = 0; i < nx; ++i) {
      std::vector<Vec3> col;
      for (const auto& r : R) col.push_back(r[i]);
      C[i] = approximate_curve(col, p, V, ys).control_points();
    }
    const std::size_t ny = C.front().size();
    std::vector<Vec3> P;
    for (std::size_t i = 0; i < nx; ++i) P.insert(P.end(), C[i].begin(), C[i].end());
    NurbsSurface s{p, p, U, V, nx, ny, std::move(P)};

    std::vector<bool> badx(bx.size() - 1, false), bady(by.size() - 1, false);
    bool any = false;
    for (std::size_t a = 0; a + 1 < bx.size(); ++a)
      for (std::size_t b = 0; b + 1 < by.size(); ++b)
        for (int q = 0; q <= 8 && !(badx[a] && bady[b]); ++q)
          for (int r = 0; r <= 8; ++r) {
            const double x = bx[a] + (bx[a + 1] - bx[a]) * q / 8.0, y = by[b] + (by[b + 1] - by[b]) * r / 8.0;
            const Vec3 pt = s.evaluate(x, y);
            if (std::abs(pt.z - f(x, y)) > tol || std::abs(pt.x - x) > tol || std::abs(pt.y - y) > tol) {
              badx[a] = bady[b] = any = true;
              break;
            }
          }
    if (!any) return s;
    if (bx.size() > max_spans || by.size() > max_spans)
      throw GeometryError("sag_surface: tolerance not reached with " + std::to_string(max_spans) + " spans");
    refine(bx, badx);
    refine(by, bady);
  }
}

NurbsSurface zernike_surface(const ZernikeSag& z, double half_width, double tolerance) {
  const double a = half_width > 0.0 ? half_width : z.norm_radius;
  return sag_surface([&](double x, double y) { return z.sag(x, y); }, -a, a, -a, a, 3, tolerance);
}

NurbsSurface xy_polynomial_surface(const XYPolynomialSag& p, double half_width, double tolerance) {
  const double a = half_width > 0.0 ? half_width : p.norm_radius;
  return sag_surface([&](double x, double y) { return p.sag(x, y); }, -a, a, -a, a, 3, tolerance);
}

}  // namespace nf
