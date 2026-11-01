// Nurbforge benchmarks: wall-clock times of representative kernel operations.
//
//   nf_bench            full run (best of several repetitions per case)
//   nf_bench --quick    one repetition each (smoke test in ctest)
//
// Prints a Markdown table. Build with -DNF_BUILD_BENCHMARKS=ON (Release).

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <functional>
#include <numbers>
#include <string>
#include <vector>

#include "nurbforge/nurbforge.hpp"

using namespace nf;

namespace {

constexpr double pi = std::numbers::pi;

double sink = 0.0;  // results are accumulated here so that the work is not optimised away

Shell sphere_shell(const Vec3& c, double r) {
  return Shell::sew({TrimmedSurface(revolve(make_arc(c, {0, 0, -1}, {1, 0, 0}, r, 0.0, pi), c, {0, 0, 1}, 2 * pi))});
}

Shell box_shell(const Vec3& lo, const Vec3& hi) {
  const Vec3 d = hi - lo;
  return Shell::sew({TrimmedSurface(make_rectangle(lo, {d.x, 0, 0}, {0, d.y, 0})),
                     TrimmedSurface(make_rectangle({lo.x, lo.y, hi.z}, {d.x, 0, 0}, {0, d.y, 0})),
                     TrimmedSurface(make_rectangle(lo, {d.x, 0, 0}, {0, 0, d.z})),
                     TrimmedSurface(make_rectangle({lo.x, hi.y, lo.z}, {d.x, 0, 0}, {0, 0, d.z})),
                     TrimmedSurface(make_rectangle(lo, {0, d.y, 0}, {0, 0, d.z})),
                     TrimmedSurface(make_rectangle({hi.x, lo.y, lo.z}, {0, d.y, 0}, {0, 0, d.z}))});
}

Shell cylinder_shell(double r, double z0, double z1) {
  const auto side = extrude(make_circle({0, 0, z0}, {0, 0, 1}, r), {0, 0, z1 - z0});
  const auto lid = [&](double z) {
    return TrimmedSurface(make_rectangle({-r, -r, z}, {2 * r, 0, 0}, {0, 2 * r, 0}), circle_loop(0.5, 0.5, 0.5));
  };
  return Shell::sew({TrimmedSurface(side), lid(z0), lid(z1)});
}

struct Case {
  std::string name;
  std::string unit;  // what one repetition does
  std::function<void()> run;
};

}  // namespace

int main(int argc, char** argv) {
  const bool quick = argc > 1 && std::strcmp(argv[1], "--quick") == 0;
  const int reps = quick ? 1 : 5;

  // fixtures
  std::vector<Vec3> pts;
  for (int i = 0; i < 1000; ++i) {
    const double t = 8.0 * pi * i / 999.0;
    pts.push_back({std::cos(t), std::sin(t), 0.1 * t});
  }
  const NurbsCurve helix = interpolate_curve(pts, 3);
  const NurbsSurface sphere = revolve(make_arc({0, 0, 0}, {0, 0, -1}, {1, 0, 0}, 1.0, 0.0, pi), {0, 0, 0}, {0, 0, 1}, 2 * pi);
  const NurbsSurface sphere2 =
      revolve(make_arc({1, 0, 0}, {0, 0, -1}, {1, 0, 0}, 1.0, 0.0, pi), {1, 0, 0}, {0, 0, 1}, 2 * pi);
  const Shell s1 = sphere_shell({0, 0, 0}, 1.0), s2 = sphere_shell({1, 0, 0}, 1.0);
  const Shell cube = box_shell({-1, -1, -1}, {1, 1, 1}), drill = cylinder_shell(0.3, -2, 2);
  const Shell box_a = box_shell({0, 0, 0}, {1, 1, 1}), box_b = box_shell({0.5, 0, 0}, {1.5, 1, 1});
  const Shell lens = boolean(s1, s2, BooleanOp::Intersection);

  const std::vector<Case> cases = {
      {"curve evaluation (cubic, 1000 ctrl)", "1e5 points",
       [&] {
         const double a = helix.domain_start(), b = helix.domain_end();
         for (int i = 0; i < 100000; ++i) sink += helix.evaluate(a + (b - a) * i / 99999.0).z;
       }},
      {"surface evaluation (rational sphere)", "1e5 points",
       [&] {
         for (int i = 0; i < 100000; ++i) sink += sphere.evaluate(std::fmod(i * 0.618, 1.0), (i % 997) / 997.0).x;
       }},
      {"closest point on sphere", "1000 points",
       [&] {
         for (int i = 0; i < 1000; ++i) sink += closest_point(sphere, {2 * std::cos(i * 1.0), 2 * std::sin(i * 1.0), 0.3}).u;
       }},
      {"curve interpolation (1000 points)", "1 fit", [&] { sink += static_cast<double>(interpolate_curve(pts, 3).control_points().size()); }},
      {"curve approximation to 1e-4", "1 fit",
       [&] { sink += static_cast<double>(approximate_curve_tol(pts, 3, 1e-4).control_points().size()); }},
      {"surface/surface intersection (two spheres)", "1 intersection",
       [&] { sink += static_cast<double>(intersect(sphere, sphere2).curves.size()); }},
      {"sewing a cube", "1 shell", [&] { sink += static_cast<double>(box_shell({0, 0, 0}, {1, 1, 1}).faces().size()); }},
      {"point in shell (sphere)", "100 points",
       [&] {
         for (int i = 0; i < 100; ++i) sink += point_in_shell(s1, {0.5 * std::cos(i * 1.0), 0.3, 0.5 * std::sin(i * 1.0)}) ? 1.0 : 0.0;
       }},
      {"Boolean intersection (two spheres)", "1 Boolean",
       [&] { sink += static_cast<double>(boolean(s1, s2, BooleanOp::Intersection).faces().size()); }},
      {"Boolean difference (drilled cube)", "1 Boolean",
       [&] { sink += static_cast<double>(boolean(cube, drill, BooleanOp::Difference).faces().size()); }},
      {"Boolean union (coincident boxes)", "1 Boolean",
       [&] { sink += static_cast<double>(boolean(box_a, box_b, BooleanOp::Union).faces().size()); }},
      {"volume of the lens (rel. 1e-6)", "1 volume", [&] { sink += lens.volume(1e-6); }},
  };

  std::printf("| case | per | time [ms] |\n|---|---|---:|\n");
  for (const auto& c : cases) {
    double best = 1e300;
    for (int r = 0; r < reps; ++r) {
      const auto t0 = std::chrono::steady_clock::now();
      c.run();
      best = std::min(best, std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count());
    }
    std::printf("| %s | %s | %.2f |\n", c.name.c_str(), c.unit.c_str(), best);
  }
  std::printf("\n(checksum %g)\n", sink);
  return 0;
}
