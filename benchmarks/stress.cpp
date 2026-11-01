// Boolean robustness stress test across surface families.
//
//   nf_stress [cases per family] [family filter]     e.g.  nf_stress 20  |  nf_stress 5 torus
//
// Random but reproducible configurations (fixed seeds, every random number drawn in a separate
// statement, so all compilers generate the same cases). Oracle, independent of analytic volume
// formulas: every result is closed and consistently oriented, and
//     V(A u B) + V(A n B) = V(A) + V(B),   V(A - B) = V(A) - V(A n B)   (relative 1e-5).
// Exit code 1 if any case fails. ctest runs one case per family (stress_smoke).

#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <numbers>
#include <string>
#include <utility>
#include <vector>

#include "nurbforge/nurbforge.hpp"

using namespace nf;

namespace {

constexpr double pi = std::numbers::pi;

struct Rng {
  std::uint64_t s;
  double operator()(double lo, double hi) {
    s = s * 6364136223846793005ULL + 1442695040888963407ULL;
    return lo + (hi - lo) * static_cast<double>(s >> 11) / static_cast<double>(1ULL << 53);
  }
};

std::uint64_t fnv1a(const char* text) {
  std::uint64_t h = 1469598103934665603ULL;
  for (const char* c = text; *c; ++c) h = (h ^ static_cast<unsigned char>(*c)) * 1099511628211ULL;
  return h;
}

Shell sphere(const Vec3& c, double r) {
  return Shell::sew({TrimmedSurface(revolve(make_arc(c, {0, 0, -1}, {1, 0, 0}, r, 0.0, pi), c, {0, 0, 1}, 2 * pi))});
}

Shell box(const Vec3& lo, const Vec3& hi) {
  const Vec3 d = hi - lo;
  return Shell::sew({TrimmedSurface(make_rectangle(lo, {d.x, 0, 0}, {0, d.y, 0})),
                     TrimmedSurface(make_rectangle({lo.x, lo.y, hi.z}, {d.x, 0, 0}, {0, d.y, 0})),
                     TrimmedSurface(make_rectangle(lo, {d.x, 0, 0}, {0, 0, d.z})),
                     TrimmedSurface(make_rectangle({lo.x, hi.y, lo.z}, {d.x, 0, 0}, {0, 0, d.z})),
                     TrimmedSurface(make_rectangle(lo, {0, d.y, 0}, {0, 0, d.z})),
                     TrimmedSurface(make_rectangle({hi.x, lo.y, lo.z}, {0, d.y, 0}, {0, 0, d.z}))});
}

Shell cylinder(const Vec3& base, double r, double h) {  // axis z
  const auto side = extrude(make_circle(base, {0, 0, 1}, r), {0, 0, h});
  const auto lid = [&](double z) {
    return TrimmedSurface(make_rectangle({base.x - r, base.y - r, z}, {2 * r, 0, 0}, {0, 2 * r, 0}),
                          circle_loop(0.5, 0.5, 0.5));
  };
  return Shell::sew({TrimmedSurface(side), lid(base.z), lid(base.z + h)});
}

Shell torus(const Vec3& c, double R, double r) {  // axis z
  const auto profile = make_circle({c.x + R, c.y, c.z}, {0, 1, 0}, r);
  return Shell::sew({TrimmedSurface(revolve(profile, c, {0, 0, 1}, 2 * pi))});
}

Shell rotated(const Shell& s, const Vec3& axis, double angle) {
  std::vector<TrimmedSurface> faces;
  for (const auto& f : s.faces())
    faces.emplace_back(f.surface().transformed(Transform::rotation({0, 0, 0}, axis, angle)), f.outer(), f.holes());
  return Shell::sew(std::move(faces));
}

Vec3 draw_vec(Rng& g, double lo, double hi) {
  const double x = g(lo, hi);
  const double y = g(lo, hi);
  const double z = g(lo, hi);
  return {x, y, z};
}

struct Family {
  const char* name;
  std::function<std::pair<Shell, Shell>(Rng&)> make;
};

const std::vector<Family>& families() {
  static const std::vector<Family> fams = {
      {"sphere/sphere",
       [](Rng& g) {
         const double r1 = g(0.6, 1.5);
         const double r2 = g(0.6, 1.5);
         Vec3 c = draw_vec(g, -1.5, 1.5);
         const double d = norm(c);
         if (d >= r1 + r2 - 0.05 || d <= std::abs(r1 - r2) + 0.05) c = normalized(c) * (0.5 * (r1 + r2 + std::abs(r1 - r2)));
         return std::pair{sphere({0, 0, 0}, r1), sphere(c, r2)};
       }},
      {"sphere/cylinder",
       [](Rng& g) {
         const double R = g(0.8, 1.2);
         const double h = g(1.0, 3.0);
         const double r = g(0.2, 0.6);
         const double bx = g(-0.6, 0.6);
         const double by = g(-0.6, 0.6);
         const double bz = g(-2.0, -0.2);
         return std::pair{sphere({0, 0, 0}, R), cylinder({bx, by, bz}, r, h)};
       }},
      {"sphere/box",
       [](Rng& g) {
         const Vec3 lo = draw_vec(g, -1.5, 0.3);
         const Vec3 size = draw_vec(g, 0.5, 1.8);
         return std::pair{sphere({0, 0, 0}, 1.0), box(lo, lo + size)};
       }},
      {"cylinder/box",
       [](Rng& g) {
         const double lx = g(-1.0, 0.2);
         const double ly = g(-1.0, 0.2);
         const double lz = g(-0.5, 1.0);
         const Vec3 lo{lx, ly, lz};
         const Vec3 size = draw_vec(g, 0.4, 1.5);
         return std::pair{cylinder({0, 0, 0}, 0.7, 2.0), box(lo, lo + size)};
       }},
      {"torus/box",
       [](Rng& g) {
         const double lx = g(-1.8, 0.5);
         const double ly = g(-1.8, 0.5);
         const double lz = g(-0.6, 0.0);
         const double sx = g(0.5, 1.5);
         const double sy = g(0.5, 1.5);
         const double sz = g(0.3, 1.0);
         const Vec3 lo{lx, ly, lz};
         return std::pair{torus({0, 0, 0}, 1.0, 0.35), box(lo, lo + Vec3{sx, sy, sz})};
       }},
      {"torus/sphere",
       [](Rng& g) {
         const double cx = g(0.3, 1.5);
         const double cy = g(-0.5, 0.5);
         const double cz = g(-0.3, 0.3);
         const double r = g(0.3, 0.7);
         return std::pair{torus({0, 0, 0}, 1.0, 0.35), sphere({cx, cy, cz}, r)};
       }},
      {"cylinder/cylinder",
       [](Rng& g) {
         const double r = g(0.2, 0.5);
         const double bx = g(-0.3, 0.3);
         const double by = g(-0.3, 0.3);
         const Shell cross = rotated(cylinder({0, 0, -1.5}, r, 3.0), {0, 1, 0}, pi / 2);  // along x
         return std::pair{cylinder({bx, by, -1.0}, 0.6, 2.0), cross};
       }},
  };
  return fams;
}

}  // namespace

int main(int argc, char** argv) {
  const int n = argc > 1 ? std::atoi(argv[1]) : 10;
  const char* only = argc > 2 ? argv[2] : nullptr;
  int total = 0, failed = 0;
  double worst = 0.0, seconds = 0.0;
  for (const auto& fam : families()) {
    if (only && std::strstr(fam.name, only) == nullptr) continue;
    Rng g{fnv1a(fam.name)};
    int fam_failed = 0;
    for (int k = 0; k < n; ++k) {
      ++total;
      try {
        const auto [a, b] = fam.make(g);
        const double va = a.oriented().volume(1e-7), vb = b.oriented().volume(1e-7);
        const auto t0 = std::chrono::steady_clock::now();
        const Shell u = boolean(a, b, BooleanOp::Union);
        const Shell i = boolean(a, b, BooleanOp::Intersection);
        const Shell d = boolean(a, b, BooleanOp::Difference);
        seconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        bool ok = true;
        for (const Shell* s : {&u, &i, &d})
          if (!s->is_closed() || !s->is_consistently_oriented()) ok = false;
        if (!ok) {
          ++failed;
          ++fam_failed;
          std::printf("%s #%d: open or inconsistently oriented result\n", fam.name, k);
          continue;
        }
        const double vu = u.volume(1e-7), vi = i.volume(1e-7), vd = d.volume(1e-7);
        const double e1 = std::abs(vu + vi - va - vb) / (va + vb), e2 = std::abs(vd - (va - vi)) / va;
        worst = std::max({worst, e1, e2});
        if (e1 > 1e-5 || e2 > 1e-5 || vi < -1e-9) {
          ++failed;
          ++fam_failed;
          std::printf("%s #%d: volume identity error %.2e / %.2e\n", fam.name, k, e1, e2);
        }
      } catch (const std::exception& e) {
        ++failed;
        ++fam_failed;
        std::printf("%s #%d: exception: %s\n", fam.name, k, e.what());
      }
    }
    std::printf("%-18s %d/%d failed\n", fam.name, fam_failed, n);
    std::fflush(stdout);
  }
  std::printf("total %d cases (3 Booleans each), %d failed, worst identity error %.2e, %.2f s per case\n", total,
              failed, worst, total > 0 ? seconds / total : 0.0);
  return failed == 0 ? 0 : 1;
}
