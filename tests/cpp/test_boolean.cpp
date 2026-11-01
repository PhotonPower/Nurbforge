#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <numbers>
#include <vector>

#include "nurbforge/core/error.hpp"
#include "nurbforge/geometry/primitives.hpp"
#include "nurbforge/modeling/construct.hpp"
#include "nurbforge/topology/boolean.hpp"

using namespace nf;
constexpr double pi = std::numbers::pi;

namespace {

bool rel_close(double a, double b, double rel) { return std::abs(a - b) <= rel * std::abs(b); }

Shell sphere_shell(const Vec3& c, double R) {
  return Shell::sew({TrimmedSurface(revolve(make_arc(c, {0, 0, -1}, {1, 0, 0}, R, 0.0, pi), c, {0, 0, 1}, 2 * pi))});
}

Shell box_shell(double h) {  // cube [-h, h]^3
  const double s = 2 * h;
  return Shell::sew({TrimmedSurface(make_rectangle({-h, -h, -h}, {s, 0, 0}, {0, s, 0})),
                     TrimmedSurface(make_rectangle({-h, -h, h}, {s, 0, 0}, {0, s, 0})),
                     TrimmedSurface(make_rectangle({-h, -h, -h}, {s, 0, 0}, {0, 0, s})),
                     TrimmedSurface(make_rectangle({-h, h, -h}, {s, 0, 0}, {0, 0, s})),
                     TrimmedSurface(make_rectangle({-h, -h, -h}, {0, s, 0}, {0, 0, s})),
                     TrimmedSurface(make_rectangle({h, -h, -h}, {0, s, 0}, {0, 0, s}))});
}

Shell box_shell(const Vec3& lo, const Vec3& hi) {  // axis-aligned box
  const Vec3 d = hi - lo;
  return Shell::sew({TrimmedSurface(make_rectangle(lo, {d.x, 0, 0}, {0, d.y, 0})),
                     TrimmedSurface(make_rectangle({lo.x, lo.y, hi.z}, {d.x, 0, 0}, {0, d.y, 0})),
                     TrimmedSurface(make_rectangle(lo, {d.x, 0, 0}, {0, 0, d.z})),
                     TrimmedSurface(make_rectangle({lo.x, hi.y, lo.z}, {d.x, 0, 0}, {0, 0, d.z})),
                     TrimmedSurface(make_rectangle(lo, {0, d.y, 0}, {0, 0, d.z})),
                     TrimmedSurface(make_rectangle({hi.x, lo.y, lo.z}, {0, d.y, 0}, {0, 0, d.z}))});
}

// Union, intersection and difference volumes (an empty result is an empty shell).
void check_volumes(const Shell& a, const Shell& b, double vu, double vi, double vd) {
  const BooleanOp ops[3] = {BooleanOp::Union, BooleanOp::Intersection, BooleanOp::Difference};
  const double ref[3] = {vu, vi, vd};
  for (int o = 0; o < 3; ++o) {
    const Shell r = boolean(a, b, ops[o]);
    INFO("operation " << o);
    REQUIRE(r.is_closed());
    CHECK(r.is_consistently_oriented());
    if (ref[o] == 0.0) CHECK(r.faces().empty());
    else CHECK(rel_close(r.volume(1e-8), ref[o], 1e-6));
  }
}

Shell cylinder_shell(double R, double z0, double z1) {  // axis z
  const auto side = extrude(make_circle({0, 0, z0}, {0, 0, 1}, R), {0, 0, z1 - z0});
  const auto lid = [&](double z) {
    return TrimmedSurface(make_rectangle({-R, -R, z}, {2 * R, 0, 0}, {0, 2 * R, 0}), circle_loop(0.5, 0.5, 0.5));
  };
  return Shell::sew({TrimmedSurface(side), lid(z0), lid(z1)});
}

}  // namespace

TEST_CASE("point in shell", "[boolean]") {
  const Shell s = sphere_shell({0, 0, 0}, 1.0);
  CHECK(point_in_shell(s, {0.1, 0.2, -0.3}));
  CHECK_FALSE(point_in_shell(s, {0.9, 0.9, 0.0}));
  const Shell b = box_shell(1.0);
  CHECK(point_in_shell(b, {0.99, -0.5, 0.3}));
  CHECK_FALSE(point_in_shell(b, {1.01, 0.0, 0.0}));
}

TEST_CASE("Boolean operations of two spheres", "[boolean]") {
  const double R = 1.0, d = 1.0;
  const Shell a = sphere_shell({0, 0, 0}, R), b = sphere_shell({d, 0.0, 0.0}, R);
  const double lens = pi * (4 * R + d) * (2 * R - d) * (2 * R - d) / 12.0;
  const double ball = 4.0 / 3.0 * pi * R * R * R;
  const Shell i = boolean(a, b, BooleanOp::Intersection);
  CHECK(i.is_closed());
  CHECK(i.is_consistently_oriented());
  CHECK(rel_close(i.volume(1e-8), lens, 1e-6));
  const Shell u = boolean(a, b, BooleanOp::Union);
  CHECK(u.is_closed());
  CHECK(rel_close(u.volume(1e-8), 2 * ball - lens, 1e-6));
  const Shell m = boolean(a, b, BooleanOp::Difference);
  CHECK(m.is_closed());
  CHECK(rel_close(m.volume(1e-8), ball - lens, 1e-6));
}

TEST_CASE("cube minus an inner sphere leaves a void", "[boolean]") {
  const double r = 0.5;
  const Shell m = boolean(box_shell(1.0), sphere_shell({0.1, -0.1, 0.2}, r), BooleanOp::Difference);
  CHECK(m.is_closed());
  CHECK(m.is_consistently_oriented());
  CHECK(m.faces().size() == 7);
  CHECK(rel_close(m.volume(1e-8), 8.0 - 4.0 / 3.0 * pi * r * r * r, 1e-7));
}

TEST_CASE("drilling a cylinder through a cube", "[boolean]") {
  const double r = 0.3;
  const Shell m = boolean(box_shell(1.0), cylinder_shell(r, -2.0, 2.0), BooleanOp::Difference);
  CHECK(m.is_closed());
  CHECK(m.is_consistently_oriented());
  CHECK(rel_close(m.volume(1e-8), 8.0 - pi * r * r * 2.0, 1e-6));
  const Shell plug = boolean(box_shell(1.0), cylinder_shell(r, -2.0, 2.0), BooleanOp::Intersection);
  CHECK(plug.is_closed());
  CHECK(rel_close(plug.volume(1e-8), pi * r * r * 2.0, 1e-6));
}

TEST_CASE("Boolean input checks", "[boolean]") {
  Shell open = Shell::sew({TrimmedSurface(make_rectangle({0, 0, 0}, {1, 0, 0}, {0, 1, 0}))});
  CHECK_THROWS_AS(boolean(open, box_shell(1.0), BooleanOp::Union), GeometryError);
}

TEST_CASE("Booleans with a closed intersection curve away from the seams", "[boolean]") {
  // found by a stress test: the intersection circle crosses neither sphere's seam, so splitting
  // produces a closed cut curve (hole assignment of its clockwise twin cycle)
  const double r1 = 0.797, r2 = 1.277;
  const Vec3 c2{-1.013, 0.967, -1.294};
  const double d = norm(c2);
  const double h1 = (r2 - r1 + d) * (r2 + r1 - d) / (2 * d), h2 = (r1 - r2 + d) * (r1 + r2 - d) / (2 * d);
  const double lens = pi * h1 * h1 * (3 * r1 - h1) / 3 + pi * h2 * h2 * (3 * r2 - h2) / 3;
  const Shell a = sphere_shell({0, 0, 0}, r1), b = sphere_shell(c2, r2);
  const Shell i = boolean(a, b, BooleanOp::Intersection);
  REQUIRE(i.is_closed());
  CHECK(rel_close(i.volume(1e-7), lens, 1e-5));
  const Shell m = boolean(b, a, BooleanOp::Difference);
  REQUIRE(m.is_closed());
  CHECK(rel_close(m.volume(1e-7), 4.0 / 3.0 * pi * r2 * r2 * r2 - lens, 1e-5));
}

TEST_CASE("Booleans with coincident faces (boxes)", "[boolean][coincident]") {
  const Shell a = box_shell({0, 0, 0}, {1, 1, 1});
  SECTION("partially shared side planes") { check_volumes(a, box_shell({0.5, 0, 0}, {1.5, 1, 1}), 1.5, 0.5, 0.5); }
  SECTION("stacked, one common face") { check_volumes(a, box_shell({0, 0, 1}, {1, 1, 2}), 2.0, 0.0, 1.0); }
  SECTION("identical boxes") { check_volumes(a, box_shell({0, 0, 0}, {1, 1, 1}), 1.0, 1.0, 0.0); }
  SECTION("box flush with the top face") {
    check_volumes(a, box_shell({0.25, 0.25, 0.5}, {0.75, 0.75, 1}), 1.0, 0.125, 0.875);
  }
}

TEST_CASE("Booleans with coincident faces (cylinders)", "[boolean][coincident]") {
  const double v = pi * 0.25;  // R = 0.5, height 1
  SECTION("stacked, common lid") { check_volumes(cylinder_shell(0.5, 0, 1), cylinder_shell(0.5, 1, 2), 2 * v, 0, v); }
  SECTION("coaxial, overlapping") {
    check_volumes(cylinder_shell(0.5, 0, 1), cylinder_shell(0.5, 0.5, 1.5), 1.5 * v, 0.5 * v, 0.5 * v);
  }
}
