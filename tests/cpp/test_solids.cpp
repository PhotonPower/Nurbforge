#include <array>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <numbers>
#include <vector>

#include "nurbforge/core/error.hpp"
#include "nurbforge/geometry/primitives.hpp"
#include "nurbforge/geometry/transform.hpp"
#include "nurbforge/modeling/construct.hpp"
#include "nurbforge/topology/boolean.hpp"
#include "nurbforge/topology/solids.hpp"

using namespace nf;
using Catch::Approx;
constexpr double pi = std::numbers::pi;

namespace {

double I(const MassProperties& mp, int i, int j) {
  return mp.inertia[static_cast<std::size_t>(3 * i + j)];
}

/// Inertia tensor of a body with axial symmetry about unit axis a: Ia along a, It across.
void check_axial(const MassProperties& mp, const Vec3& a, double Ia, double It, double eps) {
  const std::array<double, 3> av{a.x, a.y, a.z};
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) {
      const double expected = It * (i == j ? 1.0 : 0.0) + (Ia - It) *
                                                              av[static_cast<std::size_t>(i)] *
                                                              av[static_cast<std::size_t>(j)];
      CHECK(I(mp, i, j) == Approx(expected).margin(eps * std::max(Ia, It)));
    }
}

}  // namespace

TEST_CASE("solids: box", "[solids]") {
  const auto b = make_box({1, 2, 3}, {2, 3, 4});
  REQUIRE(b.is_closed());
  CHECK(b.faces().size() == 6);
  const auto mp = b.mass_properties();
  CHECK(mp.volume == Approx(24.0).epsilon(1e-12));
  CHECK(mp.centroid.x == Approx(2.0).epsilon(1e-12));
  CHECK(mp.centroid.y == Approx(3.5).epsilon(1e-12));
  CHECK(mp.centroid.z == Approx(5.0).epsilon(1e-12));
  // V/12 (b^2 + c^2, a^2 + c^2, a^2 + b^2)
  CHECK(I(mp, 0, 0) == Approx(50.0).epsilon(1e-10));
  CHECK(I(mp, 1, 1) == Approx(40.0).epsilon(1e-10));
  CHECK(I(mp, 2, 2) == Approx(26.0).epsilon(1e-10));
  CHECK(std::abs(I(mp, 0, 1)) < 1e-9);
  CHECK(mp.principal_moments[0] == Approx(26.0).epsilon(1e-10));
  CHECK(std::abs(mp.principal_axes[0].z) == Approx(1.0).epsilon(1e-10));
  CHECK_THROWS_AS(make_box({0, 0, 0}, {1, 0, 1}), GeometryError);
}

TEST_CASE("solids: sphere, cylinder, cone, torus", "[solids]") {
  SECTION("sphere") {
    const double R = 2;
    const auto mp = make_sphere({1, -1, 0.5}, R).mass_properties();
    const double V = 4.0 / 3.0 * pi * R * R * R;
    CHECK(mp.volume == Approx(V).epsilon(1e-9));
    CHECK(distance(mp.centroid, Vec3{1, -1, 0.5}) < 1e-8);
    check_axial(mp, {0, 0, 1}, 0.4 * V * R * R, 0.4 * V * R * R, 1e-8);
  }
  SECTION("cylinder on an oblique axis") {
    const double r = 1.5, h = 4;
    const Vec3 a = normalized(Vec3{1, 1, 1});
    const auto c = make_cylinder({0.3, 0, -1}, {1, 1, 1}, r, h);
    REQUIRE(c.is_closed());
    const auto mp = c.mass_properties();
    const double V = pi * r * r * h;
    CHECK(mp.volume == Approx(V).epsilon(1e-9));
    CHECK(distance(mp.centroid, Vec3{0.3, 0, -1} + 0.5 * h * a) < 1e-8);
    check_axial(mp, a, 0.5 * V * r * r, V * (3 * r * r + h * h) / 12, 1e-8);
    // the distinct principal moment belongs to the axis
    CHECK(std::abs(dot(mp.principal_axes[0], a)) == Approx(1.0).epsilon(1e-8));
  }
  SECTION("cone with apex and frustum") {
    const double r = 2, h = 3;
    const auto cone = make_cone({0, 0, 0}, {0, 0, 1}, r, 0, h);
    REQUIRE(cone.is_closed());
    const auto mp = cone.mass_properties();
    const double V = pi * r * r * h / 3;
    CHECK(mp.volume == Approx(V).epsilon(1e-9));
    CHECK(mp.centroid.z == Approx(h / 4).epsilon(1e-9));
    check_axial(mp, {0, 0, 1}, 0.3 * V * r * r, V * (3 * r * r / 20 + 3 * h * h / 80), 1e-8);
    const double r1 = 1;
    const auto fr = make_cone({0, 0, 0}, {0, 0, 1}, r, r1, h).mass_properties();
    const double s = r * r + r * r1 + r1 * r1;
    CHECK(fr.volume == Approx(pi * h * s / 3).epsilon(1e-9));
    CHECK(fr.centroid.z == Approx(h * (r * r + 2 * r * r1 + 3 * r1 * r1) / (4 * s)).epsilon(1e-9));
    CHECK_THROWS_AS(make_cone({0, 0, 0}, {0, 0, 1}, 0, 0, 1), GeometryError);
  }
  SECTION("torus") {
    const double R = 3, r = 1;
    const Vec3 a = normalized(Vec3{0, 1, 2});
    const auto mp = make_torus({1, 2, 3}, a, R, r).mass_properties();
    const double V = 2 * pi * pi * R * r * r;
    CHECK(mp.volume == Approx(V).epsilon(1e-9));
    CHECK(distance(mp.centroid, Vec3{1, 2, 3}) < 1e-8);
    check_axial(mp, a, V * (R * R + 0.75 * r * r), V * (0.5 * R * R + 0.625 * r * r), 1e-8);
    CHECK_THROWS_AS(make_torus({0, 0, 0}, {0, 0, 1}, 1, 1), GeometryError);
  }
}

TEST_CASE("solids: hollow sphere and Boolean results", "[solids]") {
  SECTION("void subtracts") {
    const double R = 2, r = 1;
    const auto outer = make_sphere({0, 0, 0}, R), inner = make_sphere({0, 0, 0}, r);
    std::vector<TrimmedSurface> faces = outer.faces();
    for (const auto& f : inner.faces()) faces.push_back(f);
    const auto mp = Shell::sew(faces).oriented().mass_properties();
    CHECK(mp.volume == Approx(4.0 / 3.0 * pi * (R * R * R - r * r * r)).epsilon(1e-9));
    const double Ia = 0.4 * 4.0 / 3.0 * pi * (std::pow(R, 5) - std::pow(r, 5));
    check_axial(mp, {0, 0, 1}, Ia, Ia, 1e-8);
  }
  SECTION("box with a through hole") {
    const auto part = boolean(make_box({-2, -2, 0}, {4, 4, 1}),
                              make_cylinder({0, 0, -1}, {0, 0, 1}, 1, 3), BooleanOp::Difference);
    const auto mp = part.mass_properties();
    CHECK(mp.volume == Approx(16 - pi).epsilon(1e-7));
    CHECK(mp.centroid.z == Approx(0.5).epsilon(1e-7));
    // I_zz = box (V/12 (a^2 + b^2)) minus cylinder (pi r^2 h r^2 / 2)
    CHECK(I(mp, 2, 2) == Approx(16.0 / 12 * 32 - pi / 2).epsilon(1e-7));
  }
}

TEST_CASE("solids: transformed shells", "[solids]") {
  const auto b = make_box({0, 0, 0}, {1, 2, 3});
  const auto r = Transform::translation({5, 0, 0}) * Transform::rotation({0, 0, 0}, {1, 1, 0}, 0.7);
  const auto mb = b.mass_properties();
  const auto mr = b.transformed(r).mass_properties();
  CHECK(mr.volume == Approx(6.0).epsilon(1e-10));
  CHECK(distance(mr.centroid, r.apply(mb.centroid)) < 1e-9);
  // I' = R I R^T
  const auto& M = r.matrix();
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) {
      double e = 0;
      for (int k = 0; k < 3; ++k)
        for (int l = 0; l < 3; ++l)
          e += M[static_cast<std::size_t>(4 * i + k)] * I(mb, k, l) *
               M[static_cast<std::size_t>(4 * j + l)];
      CHECK(I(mr, i, j) == Approx(e).margin(1e-9));
    }
  // mirrors keep the normals outward, scaling scales the volume
  const auto mirrored = b.transformed(Transform::mirror({0, 0, 0}, {1, 0, 0}));
  CHECK(mirrored.volume() == Approx(6.0).epsilon(1e-10));
  CHECK(mirrored.mass_properties().centroid.x == Approx(-0.5).epsilon(1e-10));
  CHECK(make_sphere({0, 0, 0}, 1).transformed(Transform::scaling({0, 0, 0}, 2.0)).volume() ==
        Approx(32.0 / 3.0 * pi).epsilon(1e-9));
}

TEST_CASE("solids: mass properties need a closed outward shell", "[solids]") {
  const auto open = Shell::sew({TrimmedSurface(make_rectangle({0, 0, 0}, {1, 0, 0}, {0, 1, 0}))});
  CHECK_THROWS_AS(open.mass_properties(), GeometryError);
  const auto box = make_box({0, 0, 0}, {1, 1, 1});
  std::vector<TrimmedSurface> inward;
  for (const auto& f : box.faces()) inward.push_back(f.flipped());
  CHECK_THROWS_AS(Shell::sew(inward).mass_properties(), GeometryError);
}
