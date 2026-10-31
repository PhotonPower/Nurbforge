// Builds a spherical lens cap by revolving a circular arc and prints a few points.
#include <cmath>
#include <iostream>
#include <numbers>

#include "nurbforge/nurbforge.hpp"

int main() {
  using namespace nf;
  const double R = 50.0;               // radius of curvature [mm]
  const double half_aperture = 12.5;   // [mm]
  const double alpha = std::asin(half_aperture / R);

  // Meridional arc in the x-z plane, starting on the optical axis (z).
  const auto meridian = make_arc({0, 0, -R}, {0, 0, 1}, {1, 0, 0}, R, 0.0, alpha);
  const auto cap = revolve(meridian, {0, 0, 0}, {0, 0, 1}, 2 * std::numbers::pi);

  std::cout << "Nurbforge " << NURBFORGE_VERSION << " - spherical cap\n";
  std::cout << "control net: " << cap.num_u() << " x " << cap.num_v() << '\n';
  for (double v : {0.0, 0.5, 1.0}) {
    const auto p = cap.evaluate(0.125, v);
    std::cout << "S(0.125, " << v << ") = " << p << "  |p - c| = " << distance(p, {0, 0, -R}) << '\n';
  }
}
