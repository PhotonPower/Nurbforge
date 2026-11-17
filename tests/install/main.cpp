// Smoke test of the installed package: build a sphere and check its volume.
#include <cmath>
#include <cstdio>
#include <numbers>

#include "nurbforge/topology/solids.hpp"

int main() {
  const nf::Shell ball = nf::make_sphere({0.0, 0.0, 0.0}, 1.0);
  const double v = ball.volume();
  std::printf("sphere volume %.12f\n", v);
  return std::abs(v - 4.0 / 3.0 * std::numbers::pi) < 1e-6 ? 0 : 1;
}
