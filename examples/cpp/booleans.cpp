// Closed shells and Boolean operations: two overlapping spheres and a drilled cube, with the
// volumes compared against the analytic values.
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <string>

#include "nurbforge/nurbforge.hpp"

namespace {

using namespace nf;
constexpr double pi = std::numbers::pi;

// Sphere from a revolved semicircle: a single face, closed by its seam edge and two poles.
Shell sphere_shell(const Vec3& c, double r) {
  return Shell::sew({TrimmedSurface(revolve(make_arc(c, {0, 0, -1}, {1, 0, 0}, r, 0.0, pi), c, {0, 0, 1}, 2 * pi))});
}

// Cube [-h, h]^3 from six rectangles; oriented() makes all normals point outwards.
Shell box_shell(double h) {
  const double s = 2 * h;
  return Shell::sew({TrimmedSurface(make_rectangle({-h, -h, -h}, {s, 0, 0}, {0, s, 0})),
                     TrimmedSurface(make_rectangle({-h, -h, h}, {s, 0, 0}, {0, s, 0})),
                     TrimmedSurface(make_rectangle({-h, -h, -h}, {s, 0, 0}, {0, 0, s})),
                     TrimmedSurface(make_rectangle({-h, h, -h}, {s, 0, 0}, {0, 0, s})),
                     TrimmedSurface(make_rectangle({-h, -h, -h}, {0, s, 0}, {0, 0, s})),
                     TrimmedSurface(make_rectangle({h, -h, -h}, {0, s, 0}, {0, 0, s}))})
      .oriented();
}

// Cylinder along z: extruded circle plus two square lids trimmed to the inscribed circle.
Shell cylinder_shell(double r, double z0, double z1) {
  const auto side = extrude(make_circle({0, 0, z0}, {0, 0, 1}, r), {0, 0, z1 - z0});
  const auto lid = [&](double z) {
    return TrimmedSurface(make_rectangle({-r, -r, z}, {2 * r, 0, 0}, {0, 2 * r, 0}), circle_loop(0.5, 0.5, 0.5));
  };
  return Shell::sew({TrimmedSurface(side), lid(z0), lid(z1)}).oriented();
}

void report(const std::string& name, const Shell& s, double expected) {
  const double v = s.volume(1e-7);
  std::cout << std::left << std::setw(22) << name << std::right << s.faces().size() << " faces, "
            << (s.is_closed() ? "closed" : "OPEN") << ", volume " << std::setprecision(10) << v
            << " (analytic " << expected << ", rel. error " << std::setprecision(2)
            << std::abs(v - expected) / expected << ")\n";
}

}  // namespace

int main() {
  const double R = 1.0, d = 1.0;
  const Shell a = sphere_shell({0, 0, 0}, R), b = sphere_shell({d, 0, 0}, R);
  const double ball = 4.0 / 3.0 * pi * R * R * R;
  const double lens = pi * (4 * R + d) * (2 * R - d) * (2 * R - d) / 12.0;

  std::cout << "Nurbforge " << NURBFORGE_VERSION << " - Boolean operations\n";
  report("spheres: union", boolean(a, b, BooleanOp::Union), 2 * ball - lens);
  report("spheres: intersection", boolean(a, b, BooleanOp::Intersection), lens);
  report("spheres: difference", boolean(a, b, BooleanOp::Difference), ball - lens);

  const double r = 0.3;
  const Shell drilled = boolean(box_shell(1.0), cylinder_shell(r, -2.0, 2.0), BooleanOp::Difference);
  report("drilled cube", drilled, 8.0 - 2.0 * pi * r * r);
  std::cout << std::boolalpha << "origin inside the drilled cube: " << point_in_shell(drilled, {0, 0, 0})
            << ", (0.5, 0.5, 0) inside: " << point_in_shell(drilled, {0.5, 0.5, 0}) << '\n';
}
