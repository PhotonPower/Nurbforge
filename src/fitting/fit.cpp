#include "nurbforge/fitting/fit.hpp"

#include <algorithm>
#include <cmath>
#include <string>

#include "nurbforge/core/error.hpp"
#include "nurbforge/nurbs/basis.hpp"

namespace nf {

namespace {

// ------------------------------------------------------------------ dense LU (partial pivoting)

class LU {
 public:
  LU(std::vector<double> a, std::size_t n) : a_(std::move(a)), piv_(n), n_(n) {
    double amax = 0.0;
    for (double x : a_) amax = std::max(amax, std::abs(x));
    for (std::size_t k = 0; k < n_; ++k) {
      std::size_t p = k;
      for (std::size_t i = k + 1; i < n_; ++i)
        if (std::abs(at(i, k)) > std::abs(at(p, k))) p = i;
      if (std::abs(at(p, k)) <= 1e-14 * amax)
        throw GeometryError("fitting: singular system (check the data points and the degree)");
      piv_[k] = p;
      if (p != k)
        for (std::size_t j = 0; j < n_; ++j) std::swap(at(k, j), at(p, j));
      for (std::size_t i = k + 1; i < n_; ++i) {
        at(i, k) /= at(k, k);
        for (std::size_t j = k + 1; j < n_; ++j) at(i, j) -= at(i, k) * at(k, j);
      }
    }
  }

  [[nodiscard]] std::vector<Vec3> solve(std::vector<Vec3> b) const {
    for (std::size_t k = 0; k < n_; ++k) std::swap(b[k], b[piv_[k]]);
    for (std::size_t i = 1; i < n_; ++i)
      for (std::size_t j = 0; j < i; ++j) b[i] -= at(i, j) * b[j];
    for (std::size_t i = n_; i-- > 0;) {
      for (std::size_t j = i + 1; j < n_; ++j) b[i] -= at(i, j) * b[j];
      b[i] = b[i] / at(i, i);
    }
    return b;
  }

 private:
  double& at(std::size_t i, std::size_t j) { return a_[i * n_ + j]; }
  [[nodiscard]] double at(std::size_t i, std::size_t j) const { return a_[i * n_ + j]; }
  std::vector<double> a_;
  std::vector<std::size_t> piv_;
  std::size_t n_;
};

// Full row (N_0(u), ..., N_{n_ctrl-1}(u)).
std::vector<double> basis_row(std::size_t p, const KnotVector& U, std::size_t n_ctrl, double u) {
  std::vector<double> row(n_ctrl, 0.0);
  const std::size_t span = find_span(n_ctrl - 1, p, u, U);
  const auto N = basis_functions(span, u, p, U);
  for (std::size_t j = 0; j <= p; ++j) row[span - p + j] = N[j];
  return row;
}

// Raw parameters; empty if the points are all equal (degenerate row).
std::vector<double> raw_parameters(const std::vector<Vec3>& Q, Parametrization method) {
  const std::size_t m = Q.size() - 1;
  std::vector<double> t(m + 1, 0.0);
  if (method == Parametrization::Uniform) {
    for (std::size_t k = 0; k <= m; ++k) t[k] = static_cast<double>(k) / static_cast<double>(m);
    return t;
  }
  for (std::size_t k = 1; k <= m; ++k) {
    const double d = distance(Q[k], Q[k - 1]);
    t[k] = t[k - 1] + (method == Parametrization::Centripetal ? std::sqrt(d) : d);
  }
  if (t[m] <= 0.0) return {};
  for (double& x : t) x /= t[m];
  t[m] = 1.0;
  return t;
}

// Interpolation knots by averaging (Eq. 9.8).
KnotVector interpolation_knots(const std::vector<double>& ub, std::size_t p) {
  const std::size_t n = ub.size() - 1;
  KnotVector U(n + p + 2, ub.front());
  for (std::size_t i = n + 1; i < U.size(); ++i) U[i] = ub.back();
  for (std::size_t j = 1; j + p <= n; ++j) {
    double s = 0.0;
    for (std::size_t i = j; i < j + p; ++i) s += ub[i];
    U[j + p] = s / static_cast<double>(p);
  }
  return U;
}

// Approximation knots (Eq. 9.68/9.69): every knot span contains data parameters.
KnotVector approximation_knots(const std::vector<double>& ub, std::size_t p, std::size_t n_ctrl) {
  const std::size_t m = ub.size() - 1, n = n_ctrl - 1;
  KnotVector U(n + p + 2, 0.0);
  for (std::size_t i = n + 1; i < U.size(); ++i) U[i] = 1.0;
  const double d = static_cast<double>(m + 1) / static_cast<double>(n - p + 1);
  for (std::size_t j = 1; j + p <= n; ++j) {
    const double jd = static_cast<double>(j) * d;
    const auto i = static_cast<std::size_t>(jd);
    const double alpha = jd - static_cast<double>(i);
    U[p + j] = (1.0 - alpha) * ub[i - 1] + alpha * ub[i];
  }
  return U;
}

// Least squares with interpolated end points (Eq. 9.63–9.67); matrix factored once.
class LeastSquares {
 public:
  LeastSquares(std::size_t p, KnotVector U, const std::vector<double>& ub, std::size_t n_ctrl)
      : U_(std::move(U)), n_ctrl_(n_ctrl) {
    const std::size_t m = ub.size() - 1, n = n_ctrl - 1;
    for (std::size_t k = 1; k < m; ++k) rows_.push_back(basis_row(p, U_, n_ctrl, ub[k]));
    if (n >= 2) {
      const std::size_t sz = n - 1;
      std::vector<double> A(sz * sz, 0.0);
      for (const auto& r : rows_)
        for (std::size_t i = 1; i < n; ++i)
          for (std::size_t j = 1; j < n; ++j) A[(i - 1) * sz + (j - 1)] += r[i] * r[j];
      lu_.emplace(std::move(A), sz);
    }
  }

  [[nodiscard]] std::vector<Vec3> fit(const std::vector<Vec3>& Q) const {
    const std::size_t m = Q.size() - 1, n = n_ctrl_ - 1;
    std::vector<Vec3> P(n_ctrl_);
    P.front() = Q.front();
    P.back() = Q.back();
    if (!lu_) return P;
    std::vector<Vec3> b(n - 1);
    for (std::size_t k = 1; k < m; ++k) {
      const auto& r = rows_[k - 1];
      const Vec3 Rk = Q[k] - r[0] * Q.front() - r[n] * Q.back();
      for (std::size_t i = 1; i < n; ++i) b[i - 1] += r[i] * Rk;
    }
    const auto X = lu_->solve(std::move(b));
    std::copy(X.begin(), X.end(), P.begin() + 1);
    return P;
  }

  [[nodiscard]] const KnotVector& knots() const { return U_; }

 private:
  KnotVector U_;
  std::size_t n_ctrl_;
  std::vector<std::vector<double>> rows_;
  std::optional<LU> lu_;
};

LU interpolation_matrix(std::size_t p, const KnotVector& U, const std::vector<double>& ub) {
  const std::size_t n = ub.size();
  std::vector<double> A;
  A.reserve(n * n);
  for (double u : ub) {
    const auto r = basis_row(p, U, n, u);
    A.insert(A.end(), r.begin(), r.end());
  }
  return {std::move(A), n};
}

// Averaged parameters over several point sequences, skipping degenerate ones.
std::vector<double> averaged_parameters(const std::vector<std::vector<Vec3>>& seqs,
                                        Parametrization method) {
  std::vector<double> sum(seqs.front().size(), 0.0);
  std::size_t used = 0;
  for (const auto& r : seqs) {
    const auto t = raw_parameters(r, method);
    if (t.empty()) continue;
    for (std::size_t k = 0; k < t.size(); ++k) sum[k] += t[k];
    ++used;
  }
  if (used == 0) throw GeometryError("fitting: all point rows are degenerate");
  for (double& x : sum) x /= static_cast<double>(used);
  sum.front() = 0.0;
  sum.back() = 1.0;
  for (std::size_t k = 1; k < sum.size(); ++k)
    if (!(sum[k] > sum[k - 1])) throw GeometryError("fitting: coincident data parameters");
  return sum;
}

void check_grid(const std::vector<Vec3>& Q, std::size_t nu, std::size_t nv, std::size_t pu,
                std::size_t pv) {
  if (Q.size() != nu * nv) throw GeometryError("fitting: grid size must be num_u * num_v");
  if (pu < 1 || pv < 1) throw GeometryError("fitting: degrees must be >= 1");
  if (nu < pu + 1 || nv < pv + 1) throw GeometryError("fitting: need at least degree+1 points per direction");
}

std::vector<std::vector<Vec3>> columns(const std::vector<Vec3>& Q, std::size_t nu, std::size_t nv) {
  std::vector<std::vector<Vec3>> c(nv, std::vector<Vec3>(nu));
  for (std::size_t i = 0; i < nu; ++i)
    for (std::size_t j = 0; j < nv; ++j) c[j][i] = Q[i * nv + j];
  return c;
}

std::vector<std::vector<Vec3>> rows(const std::vector<Vec3>& Q, std::size_t nu, std::size_t nv) {
  std::vector<std::vector<Vec3>> r(nu);
  for (std::size_t i = 0; i < nu; ++i)
    r[i].assign(Q.begin() + static_cast<std::ptrdiff_t>(i * nv),
                Q.begin() + static_cast<std::ptrdiff_t>((i + 1) * nv));
  return r;
}

}  // namespace

std::vector<double> data_parameters(const std::vector<Vec3>& points, Parametrization method) {
  if (points.size() < 2) throw GeometryError("fitting: need at least two points");
  auto t = raw_parameters(points, method);
  if (t.empty()) throw GeometryError("fitting: all points coincide");
  for (std::size_t k = 1; k < t.size(); ++k)
    if (!(t[k] > t[k - 1])) throw GeometryError("fitting: coincident consecutive points");
  return t;
}

std::vector<double> data_parameters(const std::vector<std::vector<Vec3>>& sequences,
                                    Parametrization method) {
  if (sequences.empty() || sequences.front().size() < 2)
    throw GeometryError("fitting: need at least two points per sequence");
  for (const auto& s : sequences)
    if (s.size() != sequences.front().size())
      throw GeometryError("fitting: all sequences must have the same number of points");
  return averaged_parameters(sequences, method);
}

NurbsCurve interpolate_curve(const std::vector<Vec3>& Q, std::size_t p,
                             const std::vector<double>& params) {
  if (p < 1) throw GeometryError("fitting: degree must be >= 1");
  if (params.size() != Q.size()) throw GeometryError("interpolate_curve: one parameter per point");
  if (Q.size() < p + 1) throw GeometryError("interpolate_curve: need at least degree+1 points");
  for (std::size_t k = 1; k < params.size(); ++k)
    if (!(params[k] > params[k - 1]))
      throw GeometryError("interpolate_curve: parameters must be strictly increasing");
  const KnotVector U = interpolation_knots(params, p);
  return {p, U, interpolation_matrix(p, U, params).solve(Q)};
}

NurbsCurve interpolate_curve(const std::vector<Vec3>& Q, std::size_t p, Parametrization method,
                             const std::optional<std::pair<Vec3, Vec3>>& end_derivatives) {
  if (p < 1) throw GeometryError("fitting: degree must be >= 1");
  const auto ub = data_parameters(Q, method);
  const std::size_t n = Q.size() - 1;
  if (!end_derivatives) {
    if (Q.size() < p + 1) throw GeometryError("interpolate_curve: need at least degree+1 points");
    const KnotVector U = interpolation_knots(ub, p);
    return {p, U, interpolation_matrix(p, U, ub).solve(Q)};
  }

  // with end derivatives (Eq. 9.21/9.22): n + 3 control points
  if (p < 2) throw GeometryError("interpolate_curve: end derivatives need degree >= 2");
  if (n + 3 < p + 1) throw GeometryError("interpolate_curve: too few points for this degree");
  const std::size_t nc = n + 3;
  KnotVector U(nc + p + 1, 0.0);
  for (std::size_t i = nc; i < U.size(); ++i) U[i] = 1.0;
  for (std::size_t j = 0; j + p <= n + 1; ++j) {
    double s = 0.0;
    for (std::size_t i = j; i < j + p; ++i) s += ub[i];
    U[j + p + 1] = s / static_cast<double>(p);
  }
  std::vector<double> A;
  std::vector<Vec3> b;
  const auto add_row = [&](std::vector<double> r, const Vec3& rhs) {
    A.insert(A.end(), r.begin(), r.end());
    b.push_back(rhs);
  };
  for (std::size_t k = 0; k <= n; ++k) add_row(basis_row(p, U, nc, ub[k]), Q[k]);
  std::vector<double> d0(nc, 0.0), d1(nc, 0.0);
  d0[0] = -1.0;
  d0[1] = 1.0;
  d1[nc - 2] = -1.0;
  d1[nc - 1] = 1.0;
  add_row(d0, (U[p + 1] / static_cast<double>(p)) * end_derivatives->first);
  add_row(d1, ((1.0 - U[nc - 1]) / static_cast<double>(p)) * end_derivatives->second);
  return {p, U, LU(std::move(A), nc).solve(std::move(b))};
}

NurbsCurve approximate_curve(const std::vector<Vec3>& Q, std::size_t p, std::size_t n_ctrl,
                             Parametrization method) {
  if (p < 1) throw GeometryError("fitting: degree must be >= 1");
  if (n_ctrl < p + 1 || n_ctrl > Q.size())
    throw GeometryError("approximate_curve: need degree+1 <= n_ctrl <= number of points");
  const auto ub = data_parameters(Q, method);
  const LeastSquares ls(p, approximation_knots(ub, p, n_ctrl), ub, n_ctrl);
  return {p, ls.knots(), ls.fit(Q)};
}

NurbsCurve approximate_curve(const std::vector<Vec3>& Q, std::size_t p, const KnotVector& U,
                             const std::vector<double>& params) {
  if (p < 1) throw GeometryError("fitting: degree must be >= 1");
  if (params.size() != Q.size()) throw GeometryError("approximate_curve: one parameter per point");
  const std::size_t n_ctrl = U.size() - p - 1;
  validate_knots(U, p, n_ctrl);
  if (n_ctrl > Q.size()) throw GeometryError("approximate_curve: more control points than data points");
  const LeastSquares ls(p, U, params, n_ctrl);
  return {p, U, ls.fit(Q)};
}

NurbsCurve approximate_curve_tol(const std::vector<Vec3>& Q, std::size_t p, double tolerance,
                                 Parametrization method) {
  const auto ub = data_parameters(Q, method);
  const std::size_t N = Q.size();
  if (N < p + 1) throw GeometryError("approximate_curve_tol: need at least degree+1 points");
  const auto make = [&](std::size_t n) {
    return n == N ? interpolate_curve(Q, p, method) : approximate_curve(Q, p, n, method);
  };
  const auto ok = [&](const NurbsCurve& c) {
    for (std::size_t k = 0; k < N; ++k)
      if (distance(c.evaluate(ub[k]), Q[k]) > tolerance) return false;
    return true;
  };
  // exponential search for an admissible size, then bisection for the smallest one
  std::size_t lo = p, hi = p + 1;  // lo: known to fail (or below the minimum)
  NurbsCurve best = make(hi);
  while (!ok(best)) {
    lo = hi;
    if (hi == N) return best;  // interpolation: residual zero up to rounding
    hi = std::min(2 * hi, N);
    best = make(hi);
  }
  while (hi - lo > 1) {
    const std::size_t mid = (lo + hi) / 2;
    NurbsCurve c = make(mid);
    if (ok(c)) {
      hi = mid;
      best = std::move(c);
    } else {
      lo = mid;
    }
  }
  return best;
}

NurbsSurface interpolate_surface(const std::vector<Vec3>& Q, std::size_t nu, std::size_t nv,
                                 std::size_t pu, std::size_t pv, Parametrization method) {
  check_grid(Q, nu, nv, pu, pv);
  const auto cols = columns(Q, nu, nv);
  const auto ub = averaged_parameters(cols, method);
  const auto vb = averaged_parameters(rows(Q, nu, nv), method);
  const KnotVector U = interpolation_knots(ub, pu), V = interpolation_knots(vb, pv);
  const LU Au = interpolation_matrix(pu, U, ub), Av = interpolation_matrix(pv, V, vb);
  // through the columns (along u), then through the rows of the intermediate net (along v)
  std::vector<std::vector<Vec3>> R(nv);
  for (std::size_t j = 0; j < nv; ++j) R[j] = Au.solve(cols[j]);
  std::vector<Vec3> P(nu * nv);
  for (std::size_t i = 0; i < nu; ++i) {
    std::vector<Vec3> row(nv);
    for (std::size_t j = 0; j < nv; ++j) row[j] = R[j][i];
    const auto Prow = Av.solve(std::move(row));
    std::copy(Prow.begin(), Prow.end(), P.begin() + static_cast<std::ptrdiff_t>(i * nv));
  }
  return {pu, pv, U, V, nu, nv, std::move(P)};
}

NurbsSurface approximate_surface(const std::vector<Vec3>& Q, std::size_t nu, std::size_t nv,
                                 std::size_t pu, std::size_t pv, std::size_t cu, std::size_t cv,
                                 Parametrization method) {
  check_grid(Q, nu, nv, pu, pv);
  if (cu < pu + 1 || cu > nu || cv < pv + 1 || cv > nv)
    throw GeometryError("approximate_surface: need degree+1 <= ctrl <= points per direction");
  const auto cols = columns(Q, nu, nv);
  const auto ub = averaged_parameters(cols, method);
  const auto vb = averaged_parameters(rows(Q, nu, nv), method);
  const LeastSquares lsu(pu, approximation_knots(ub, pu, cu), ub, cu);
  const LeastSquares lsv(pv, approximation_knots(vb, pv, cv), vb, cv);
  std::vector<std::vector<Vec3>> R(nv);  // R[j]: cu points
  for (std::size_t j = 0; j < nv; ++j) R[j] = lsu.fit(cols[j]);
  std::vector<Vec3> P(cu * cv);
  for (std::size_t i = 0; i < cu; ++i) {
    std::vector<Vec3> row(nv);
    for (std::size_t j = 0; j < nv; ++j) row[j] = R[j][i];
    const auto Prow = lsv.fit(row);
    std::copy(Prow.begin(), Prow.end(), P.begin() + static_cast<std::ptrdiff_t>(i * cv));
  }
  return {pu, pv, lsu.knots(), lsv.knots(), cu, cv, std::move(P)};
}

}  // namespace nf
