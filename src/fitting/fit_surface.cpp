// Surface fitting to scattered points: penalised least squares with a thin-plate energy, solved
// as a banded symmetric positive definite system, followed by parameter correction.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <optional>
#include <utility>
#include <vector>

#include "nurbforge/analysis/closest_point.hpp"
#include "nurbforge/core/error.hpp"
#include "nurbforge/fitting/fit.hpp"
#include "nurbforge/modeling/gordon.hpp"
#include "nurbforge/nurbs/basis.hpp"
#include "nurbforge/nurbs/knot_vector.hpp"

namespace nf {

namespace {

// Symmetric positive definite band matrix (lower band of width b), Cholesky factorisation.
class BandMatrix {
 public:
  BandMatrix(std::size_t n, std::size_t b) : n_(n), b_(b), a_(n * (b + 1), 0.0) {}

  [[nodiscard]] std::size_t size() const { return n_; }
  [[nodiscard]] std::size_t band() const { return b_; }
  // (i, j) with j <= i <= j + b
  double& at(std::size_t i, std::size_t j) { return a_[i * (b_ + 1) + (i - j)]; }
  [[nodiscard]] double at(std::size_t i, std::size_t j) const { return a_[i * (b_ + 1) + (i - j)]; }
  void add(std::size_t i, std::size_t j, double v) {
    if (j > i) std::swap(i, j);
    at(i, j) += v;
  }

  void factor() {
    for (std::size_t i = 0; i < n_; ++i) {
      const std::size_t j0 = i > b_ ? i - b_ : 0;
      for (std::size_t j = j0; j <= i; ++j) {
        double s = at(i, j);
        for (std::size_t k = std::max(j0, j > b_ ? j - b_ : 0); k < j; ++k)
          s -= at(i, k) * at(j, k);
        if (j == i) {
          if (!(s > 0.0)) throw GeometryError("fit_surface: singular system (increase smoothing)");
          at(i, i) = std::sqrt(s);
        } else {
          at(i, j) = s / at(j, j);
        }
      }
    }
  }

  [[nodiscard]] std::vector<Vec3> solve(std::vector<Vec3> x) const {
    for (std::size_t i = 0; i < n_; ++i) {
      Vec3 s = x[i];
      for (std::size_t k = i > b_ ? i - b_ : 0; k < i; ++k) s = s - at(i, k) * x[k];
      x[i] = s / at(i, i);
    }
    for (std::size_t i = n_; i-- > 0;) {
      Vec3 s = x[i];
      for (std::size_t k = i + 1; k < std::min(n_, i + b_ + 1); ++k) s = s - at(k, i) * x[k];
      x[i] = s / at(i, i);
    }
    return x;
  }

 private:
  std::size_t n_, b_;
  std::vector<double> a_;
};

// Gauss-Legendre nodes and weights on [-1, 1] (Newton on P_n).
std::pair<std::vector<double>, std::vector<double>> gauss_legendre(std::size_t n) {
  std::vector<double> x(n), w(n);
  for (std::size_t i = 0; i < n; ++i) {
    double z = std::cos(std::numbers::pi * (static_cast<double>(i) + 0.75) /
                        (static_cast<double>(n) + 0.5));
    double dp = 1.0;
    for (int it = 0; it < 100; ++it) {
      double p0 = 1.0, p1 = z;
      for (std::size_t k = 2; k <= n; ++k) {
        const auto kd = static_cast<double>(k);
        const double p2 = ((2.0 * kd - 1.0) * z * p1 - (kd - 1.0) * p0) / kd;
        p0 = p1;
        p1 = p2;
      }
      dp = static_cast<double>(n) * (z * p1 - p0) / (z * z - 1.0);
      const double dz = p1 / dp;
      z -= dz;
      if (std::abs(dz) < 1e-16) break;
    }
    x[i] = z;
    w[i] = 2.0 / ((1.0 - z * z) * dp * dp);
  }
  return {x, w};
}

// Gram matrices G[d](i, k) = int N_i^(d) N_k^(d) for d = 0, 1, 2 (dense n x n, banded).
std::array<std::vector<double>, 3> gram(std::size_t p, const KnotVector& U, std::size_t n) {
  std::array<std::vector<double>, 3> G;
  for (auto& g : G) g.assign(n * n, 0.0);
  const auto [gx, gw] = gauss_legendre(p + 1);
  for (std::size_t s = p; s < n; ++s) {
    const double a = U[s], b = U[s + 1];
    if (!(b > a)) continue;
    for (std::size_t q = 0; q < gx.size(); ++q) {
      const double u = 0.5 * (a + b) + 0.5 * (b - a) * gx[q];
      const double w = 0.5 * (b - a) * gw[q];
      const auto d = basis_function_derivatives(s, u, p, 2, U);
      for (std::size_t r = 0; r <= 2; ++r)
        for (std::size_t i = 0; i <= p; ++i)
          for (std::size_t k = 0; k <= p; ++k)
            G[r][(s - p + i) * n + (s - p + k)] += w * d[r][i] * d[r][k];
    }
  }
  return G;
}

struct Layout {
  std::size_t pu, pv, nu, nv;
  KnotVector U, V;
};

// Penalised least squares for the control points; fixed[i] (if set) are kept.
std::vector<Vec3> solve_fit(const Layout& L, const std::vector<Vec3>& pts,
                            const std::vector<std::array<double, 2>>& uv, double lambda,
                            const std::vector<std::optional<Vec3>>& fixed) {
  const std::size_t n = L.nu * L.nv;
  BandMatrix M(n, L.pu * L.nv + L.pv);
  std::vector<Vec3> rhs(n, Vec3{0, 0, 0});
  const double inv = 1.0 / static_cast<double>(pts.size());
  std::vector<std::size_t> idx((L.pu + 1) * (L.pv + 1));
  std::vector<double> w(idx.size());
  for (std::size_t k = 0; k < pts.size(); ++k) {
    const std::size_t su = find_span(L.nu - 1, L.pu, uv[k][0], L.U);
    const std::size_t sv = find_span(L.nv - 1, L.pv, uv[k][1], L.V);
    const auto Nu = basis_functions(su, uv[k][0], L.pu, L.U);
    const auto Nv = basis_functions(sv, uv[k][1], L.pv, L.V);
    std::size_t m = 0;
    for (std::size_t a = 0; a <= L.pu; ++a)
      for (std::size_t b = 0; b <= L.pv; ++b, ++m) {
        idx[m] = (su - L.pu + a) * L.nv + (sv - L.pv + b);
        w[m] = Nu[a] * Nv[b];
      }
    for (std::size_t r = 0; r < m; ++r) {
      rhs[idx[r]] = rhs[idx[r]] + (inv * w[r]) * pts[k];
      for (std::size_t c = 0; c <= r; ++c) M.add(idx[r], idx[c], inv * w[r] * w[c]);
    }
  }
  // thin-plate energy, scaled so that lambda weighs it relative to the data term (ratio of the
  // traces): independent of the number of points and control points
  const auto Gu = gram(L.pu, L.U, L.nu), Gv = gram(L.pv, L.V, L.nv);
  const auto energy = [&](std::size_t i, std::size_t j, std::size_t k, std::size_t l) {
    return Gu[2][i * L.nu + k] * Gv[0][j * L.nv + l] +
           2.0 * Gu[1][i * L.nu + k] * Gv[1][j * L.nv + l] +
           Gu[0][i * L.nu + k] * Gv[2][j * L.nv + l];
  };
  double tr_data = 0.0, tr_energy = 0.0;
  for (std::size_t i = 0; i < L.nu; ++i)
    for (std::size_t j = 0; j < L.nv; ++j) {
      tr_data += M.at(i * L.nv + j, i * L.nv + j);
      tr_energy += energy(i, j, i, j);
    }
  const double scale = lambda * tr_data / tr_energy;
  for (std::size_t i = 0; i < L.nu; ++i)
    for (std::size_t j = 0; j < L.nv; ++j) {
      const std::size_t I = i * L.nv + j;
      for (std::size_t k = i >= L.pu ? i - L.pu : 0; k <= i; ++k)
        for (std::size_t l = j >= L.pv ? j - L.pv : 0; l < std::min(L.nv, j + L.pv + 1); ++l) {
          const std::size_t J = k * L.nv + l;
          if (J <= I) M.at(I, J) += scale * energy(i, j, k, l);
        }
    }
  // fixed control points: move their columns to the right-hand side
  const std::size_t b = M.band();
  for (std::size_t I = 0; I < n; ++I) {
    if (!fixed[I]) continue;
    for (std::size_t R = I > b ? I - b : 0; R < std::min(n, I + b + 1); ++R) {
      if (R == I) continue;
      double& m = R > I ? M.at(R, I) : M.at(I, R);
      if (!fixed[R]) rhs[R] = rhs[R] - m * *fixed[I];
      m = 0.0;
    }
    M.at(I, I) = 1.0;
    rhs[I] = *fixed[I];
  }
  M.factor();
  return M.solve(std::move(rhs));
}

std::size_t auto_count(std::size_t requested, std::size_t p, std::size_t points) {
  if (requested > 0) return std::max(requested, p + 1);
  const auto n = static_cast<std::size_t>(std::lround(std::sqrt(static_cast<double>(points) / 4)));
  return std::clamp<std::size_t>(n, p + 1, 40);
}

double clamp01(double x) { return std::clamp(x, 0.0, 1.0); }

// Parameters on a base surface (closest points), mapped to [0, 1]^2.
std::vector<std::array<double, 2>> base_parameters(const NurbsSurface& base,
                                                   const std::vector<Vec3>& pts) {
  const auto [u0, u1] = base.domain_u();
  const auto [v0, v1] = base.domain_v();
  std::vector<std::array<double, 2>> uv;
  uv.reserve(pts.size());
  for (const auto& p : pts) {
    const SurfacePoint f = closest_point(base, p);
    uv.push_back({clamp01((f.u - u0) / (u1 - u0)), clamp01((f.v - v0) / (v1 - v0))});
  }
  return uv;
}

// Height-field parameters: coordinates in the best-fit plane over the bounding rectangle,
// aligned with the principal direction of the points.
std::vector<std::array<double, 2>> plane_parameters(const std::vector<Vec3>& pts) {
  const PlaneFit pf = fit_plane(pts);
  const Vec3 n = pf.plane.normal;
  Vec3 e1 = fit_line(pts).direction;
  e1 = e1 - dot(e1, n) * n;
  if (!(norm(e1) > 1e-9))
    e1 = std::abs(n.x) < 0.9 ? cross(n, Vec3{1, 0, 0}) : cross(n, Vec3{0, 1, 0});
  e1 = normalized(e1);
  const Vec3 e2 = cross(n, e1);
  double a0 = 1e300, a1 = -1e300, b0 = 1e300, b1 = -1e300;
  std::vector<std::array<double, 2>> uv;
  uv.reserve(pts.size());
  for (const auto& p : pts) {
    const double a = dot(p - pf.plane.point, e1), b = dot(p - pf.plane.point, e2);
    a0 = std::min(a0, a), a1 = std::max(a1, a), b0 = std::min(b0, b), b1 = std::max(b1, b);
    uv.push_back({a, b});
  }
  if (!(a1 > a0) || !(b1 > b0))
    throw GeometryError("fit_surface: the points do not span an area in their best-fit plane");
  for (auto& q : uv) q = {(q[0] - a0) / (a1 - a0), (q[1] - b0) / (b1 - b0)};
  return uv;
}

}  // namespace

FittedSurface fit_surface(const std::vector<Vec3>& points, const SurfaceFitOptions& o) {
  if (points.size() < 3) throw GeometryError("fit_surface: needs at least 3 points");
  if (o.degree_u < 1 || o.degree_v < 1)
    throw GeometryError("fit_surface: degrees must be at least 1");
  if (!(o.smoothing >= 0.0)) throw GeometryError("fit_surface: smoothing must be >= 0");
  if (!o.boundary.empty() && o.boundary.size() != 4)
    throw GeometryError("fit_surface: boundary needs four curves {c0, c1, d0, d1}");

  Layout L;
  std::vector<std::array<double, 2>> uv;
  std::vector<std::optional<Vec3>> fixed;
  if (!o.boundary.empty()) {
    NurbsSurface base = coons(o.boundary[0], o.boundary[1], o.boundary[2], o.boundary[3]);
    uv = base_parameters(base, points);  // before refinement: fewer patches to search
    if (base.degree_u() < o.degree_u) base = base.elevate_degree_u(o.degree_u - base.degree_u());
    if (base.degree_v() < o.degree_v) base = base.elevate_degree_v(o.degree_v - base.degree_v());
    const std::size_t want_u = auto_count(o.num_u, base.degree_u(), points.size());
    const std::size_t want_v = auto_count(o.num_v, base.degree_v(), points.size());
    const auto refine = [](const KnotVector& K, std::size_t have, std::size_t want) {
      std::vector<double> vals = unique_knots(K).first, X;
      while (have + X.size() < want) {
        std::size_t best = 0;
        for (std::size_t s = 1; s + 1 < vals.size(); ++s)
          if (vals[s + 1] - vals[s] > vals[best + 1] - vals[best]) best = s;
        const double mid = 0.5 * (vals[best] + vals[best + 1]);
        X.push_back(mid);
        vals.insert(vals.begin() + static_cast<std::ptrdiff_t>(best) + 1, mid);
      }
      std::sort(X.begin(), X.end());
      return X;
    };
    const auto xu = refine(base.knots_u(), base.num_u(), want_u);
    if (!xu.empty()) base = base.refine_knots_u(xu);
    const auto xv = refine(base.knots_v(), base.num_v(), want_v);
    if (!xv.empty()) base = base.refine_knots_v(xv);
    if (std::any_of(base.weights().begin(), base.weights().end(),
                    [](double w) { return w != 1.0; }))
      throw GeometryError("fit_surface: boundary curves must be non-rational");
    L = {base.degree_u(), base.degree_v(), base.num_u(),
         base.num_v(),    base.knots_u(),  base.knots_v()};
    fixed.assign(L.nu * L.nv, std::nullopt);
    for (std::size_t i = 0; i < L.nu; ++i)
      for (std::size_t j = 0; j < L.nv; ++j)
        if (i == 0 || j == 0 || i + 1 == L.nu || j + 1 == L.nv)
          fixed[i * L.nv + j] = base.control_point(i, j);
  } else {
    L.pu = o.degree_u;
    L.pv = o.degree_v;
    L.nu = auto_count(o.num_u, L.pu, points.size());
    L.nv = auto_count(o.num_v, L.pv, points.size());
    L.U = clamped_uniform_knots(L.pu, L.nu);
    L.V = clamped_uniform_knots(L.pv, L.nv);
    fixed.assign(L.nu * L.nv, std::nullopt);
    uv = o.base ? base_parameters(*o.base, points) : plane_parameters(points);
  }

  const double lambda = o.smoothing > 0.0 ? o.smoothing : 1e-12;
  const auto make = [&](std::vector<Vec3> P) {
    return NurbsSurface(L.pu, L.pv, L.U, L.V, L.nu, L.nv, std::move(P));
  };
  NurbsSurface S = make(solve_fit(L, points, uv, lambda, fixed));
  for (std::size_t it = 0; it <= o.iterations; ++it) {
    for (std::size_t k = 0; k < points.size(); ++k) {
      const SurfacePoint f = closest_point_local(S, points[k], uv[k][0], uv[k][1]);
      uv[k] = {clamp01(f.u), clamp01(f.v)};
    }
    if (it < o.iterations) S = make(solve_fit(L, points, uv, lambda, fixed));
  }

  FittedSurface r{S, 0.0, 0.0, uv};
  double sum = 0.0;
  for (std::size_t k = 0; k < points.size(); ++k) {
    const double d = distance(S.evaluate(uv[k][0], uv[k][1]), points[k]);
    r.max_error = std::max(r.max_error, d);
    sum += d * d;
  }
  r.rms_error = std::sqrt(sum / static_cast<double>(points.size()));
  return r;
}

}  // namespace nf
