#include "nurbforge/nurbs/basis.hpp"

#include <algorithm>
#include <utility>

namespace nf {

std::vector<double> basis_functions(std::size_t span, double u, std::size_t p,
                                    const KnotVector& U) {
  std::vector<double> N(p + 1, 0.0), left(p + 1, 0.0), right(p + 1, 0.0);
  N[0] = 1.0;
  for (std::size_t j = 1; j <= p; ++j) {
    left[j] = u - U[span + 1 - j];
    right[j] = U[span + j] - u;
    double saved = 0.0;
    for (std::size_t r = 0; r < j; ++r) {
      const double temp = N[r] / (right[r + 1] + left[j - r]);
      N[r] = saved + right[r + 1] * temp;
      saved = left[j - r] * temp;
    }
    N[j] = saved;
  }
  return N;
}

std::vector<std::vector<double>> basis_function_derivatives(std::size_t span, double u,
                                                            std::size_t p, std::size_t n,
                                                            const KnotVector& U) {
  using Mat = std::vector<std::vector<double>>;
  using I = long;  // signed indices keep A2.3 readable
  auto at = [](Mat& m, I i, I j) -> double& {
    return m[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)];
  };

  Mat ndu(p + 1, std::vector<double>(p + 1, 0.0));
  std::vector<double> left(p + 1, 0.0), right(p + 1, 0.0);
  ndu[0][0] = 1.0;
  for (std::size_t j = 1; j <= p; ++j) {
    left[j] = u - U[span + 1 - j];
    right[j] = U[span + j] - u;
    double saved = 0.0;
    for (std::size_t r = 0; r < j; ++r) {
      ndu[j][r] = right[r + 1] + left[j - r];
      const double temp = ndu[r][j - 1] / ndu[j][r];
      ndu[r][j] = saved + right[r + 1] * temp;
      saved = left[j - r] * temp;
    }
    ndu[j][j] = saved;
  }

  Mat ders(n + 1, std::vector<double>(p + 1, 0.0));
  for (std::size_t j = 0; j <= p; ++j) ders[0][j] = ndu[j][p];

  // derivatives of order > p vanish identically
  const I ip = static_cast<I>(p), in = static_cast<I>(std::min(n, p));
  Mat a(2, std::vector<double>(p + 1, 0.0));
  for (I r = 0; r <= ip; ++r) {
    I s1 = 0, s2 = 1;
    at(a, 0, 0) = 1.0;
    for (I k = 1; k <= in; ++k) {
      double d = 0.0;
      const I rk = r - k, pk = ip - k;
      if (r >= k) {
        at(a, s2, 0) = at(a, s1, 0) / at(ndu, pk + 1, rk);
        d = at(a, s2, 0) * at(ndu, rk, pk);
      }
      const I j1 = (rk >= -1) ? 1 : -rk;
      const I j2 = (r - 1 <= pk) ? k - 1 : ip - r;
      for (I j = j1; j <= j2; ++j) {
        at(a, s2, j) = (at(a, s1, j) - at(a, s1, j - 1)) / at(ndu, pk + 1, rk + j);
        d += at(a, s2, j) * at(ndu, rk + j, pk);
      }
      if (r <= pk) {
        at(a, s2, k) = -at(a, s1, k - 1) / at(ndu, pk + 1, r);
        d += at(a, s2, k) * at(ndu, r, pk);
      }
      at(ders, k, r) = d;
      std::swap(s1, s2);
    }
  }
  double factor = static_cast<double>(p);
  for (std::size_t k = 1; k <= n; ++k) {
    for (std::size_t j = 0; j <= p; ++j) ders[k][j] *= factor;
    factor *= (p > k) ? static_cast<double>(p - k) : 0.0;
  }
  return ders;
}

}  // namespace nf
