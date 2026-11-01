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
