#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iterator>
#include <numbers>
#include <string>
#include <vector>

#include "nurbforge/core/error.hpp"
#include "nurbforge/geometry/primitives.hpp"
#include "nurbforge/mesh/shrink_wrap.hpp"
#include "nurbforge/topology/boolean.hpp"
#include "nurbforge/topology/solids.hpp"

using namespace nf;
constexpr double pi = std::numbers::pi;

namespace {

bool rel_close(double a, double b, double rel) { return std::abs(a - b) <= rel * std::abs(b); }

// Closed, oriented, of decent shape, and with the Euler characteristic of `components` closed
// surfaces with `genus` handles in total: V - E + F = 2 (components - genus).
void check_closed(const TriangleMesh& m, int components = 1, int genus = 0) {
  CHECK(m.is_watertight());
  CHECK(m.is_consistently_oriented());
  CHECK(m.min_quality() >= 0.3);
  CHECK(static_cast<long>(m.nodes.size()) - static_cast<long>(m.triangles.size() / 2) ==
        2 * (components - genus));
}

}  // namespace

TEST_CASE("TriangleMesh: measures of a tetrahedron", "[mesh]") {
  TriangleMesh m;
  m.nodes = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
  m.triangles = {{0, 2, 1}, {0, 1, 3}, {1, 2, 3}, {0, 3, 2}};
  CHECK(m.is_watertight());
  CHECK(m.is_consistently_oriented());
  CHECK(std::abs(m.volume() - 1.0 / 6.0) < 1e-15);
  CHECK(std::abs(m.area() - (1.5 + std::sqrt(3.0) / 2.0)) < 1e-14);
  const auto q = m.quality();
  CHECK(std::abs(q[2] - 1.0) < 1e-14);                   // equilateral
  CHECK(std::abs(q[0] - std::sqrt(3.0) / 2.0) < 1e-14);  // right isosceles
  CHECK(std::abs(m.min_quality() - std::sqrt(3.0) / 2.0) < 1e-14);

  SECTION("a reversed triangle breaks the orientation, a missing one the closure") {
    TriangleMesh r = m;
    std::swap(r.triangles[0][1], r.triangles[0][2]);
    CHECK(r.is_watertight());
    CHECK_FALSE(r.is_consistently_oriented());
    r = m;
    r.triangles.pop_back();
    CHECK_FALSE(r.is_watertight());
  }
  SECTION("a degenerate triangle has quality 0") {
    TriangleMesh d;
    d.nodes = {{0, 0, 0}, {1, 0, 0}, {2, 0, 0}};
    d.triangles = {{0, 1, 2}};
    CHECK(d.quality()[0] == 0.0);
  }
  SECTION("binary STL") {
    const std::string path = "nurbforge_test_mesh.stl";
    write_stl(m, path);
    std::ifstream f(path, std::ios::binary);
    const std::vector<char> data((std::istreambuf_iterator<char>(f)),
                                 std::istreambuf_iterator<char>());
    f.close();
    std::remove(path.c_str());
    REQUIRE(data.size() == 84 + 50 * 4);
    std::uint32_t count = 0;
    std::memcpy(&count, data.data() + 80, 4);
    CHECK(count == 4);
    float first[12];
    std::memcpy(first, data.data() + 84, sizeof first);
    CHECK(first[2] == -1.0f);  // normal of the triangle in the plane z = 0
    CHECK(first[6] == 0.0f);   // second corner: node 2 = (0, 1, 0)
    CHECK(first[7] == 1.0f);
    CHECK_THROWS_AS(write_stl(m, "no/such/directory/mesh.stl"), GeometryError);
  }
}

TEST_CASE("shrink_wrap: the hull of a sphere lies on the sphere", "[mesh][shrink_wrap]") {
  const std::vector<Shell> ball{make_sphere({0, 0, 0}, 1.0)};
  for (const auto& [h, rel] : {std::pair{0.2, 0.05}, std::pair{0.1, 0.015}}) {
    const TriangleMesh m = shrink_wrap(ball, h);
    check_closed(m);
    CHECK(rel_close(m.volume(), 4.0 / 3.0 * pi, rel));
    double rmin = 1e9, rmax = 0.0;
    for (const Vec3& p : m.nodes) {
      rmin = std::min(rmin, norm(p));
      rmax = std::max(rmax, norm(p));
    }
    CHECK(rmax < 1.0 + 0.05 * h);
    CHECK(rmin > 1.0 - 0.15 * h);
    // normals point outwards
    bool outwards = true;
    for (const auto& t : m.triangles) {
      const Vec3 n = cross(m.nodes[t[1]] - m.nodes[t[0]], m.nodes[t[2]] - m.nodes[t[0]]);
      outwards = outwards && dot(n, m.nodes[t[0]]) > 0.0;
    }
    CHECK(outwards);
  }

  SECTION("offset") {
    ShrinkWrapOptions opt;
    opt.offset = 0.5;
    const TriangleMesh m = shrink_wrap(ball, 0.1, opt);
    check_closed(m);
    for (const Vec3& p : m.nodes) {
      REQUIRE(norm(p) < 1.503);
      REQUIRE(norm(p) > 1.49);
    }
  }
  SECTION("without filling, the inside of the shell stays a void") {
    ShrinkWrapOptions opt;
    opt.offset = 0.2;
    opt.fill_holes = false;
    const TriangleMesh m = shrink_wrap(ball, 0.1, opt);
    check_closed(m, 2);
    CHECK(rel_close(m.volume(), 4.0 / 3.0 * pi * (1.2 * 1.2 * 1.2 - 0.8 * 0.8 * 0.8), 0.03));
  }
  SECTION("triangle quality: only the nodes move") {
    ShrinkWrapOptions raw;
    raw.min_quality = 0.0;  // as cut out of the grid: slivers
    TriangleMesh r = shrink_wrap(ball, 0.1, raw);
    CHECK(r.is_watertight());
    CHECK(r.min_quality() < 0.01);
    const TriangleMesh m = shrink_wrap(ball, 0.1);
    CHECK(m.triangles == r.triangles);
    CHECK(rel_close(m.volume(), r.volume(), 1e-3));
    const std::size_t steps = improve_quality(r, 0.3);
    CHECK(steps > 0);
    CHECK(steps <= 50);
    CHECK(r.min_quality() >= 0.3);
    CHECK(improve_quality(r, 0.3) == 0);
  }
}

TEST_CASE("shrink_wrap: boxes, overlaps, gaps and holes", "[mesh][shrink_wrap]") {
  SECTION("a box keeps its extent") {
    const TriangleMesh m = shrink_wrap(std::vector<Shell>{make_box({0, 0, 0}, {2, 1, 1})}, 0.05);
    check_closed(m);
    CHECK(rel_close(m.volume(), 2.0, 0.01));
    Vec3 lo{1e9, 1e9, 1e9}, hi{-1e9, -1e9, -1e9};
    for (const Vec3& p : m.nodes) {
      lo = {std::min(lo.x, p.x), std::min(lo.y, p.y), std::min(lo.z, p.z)};
      hi = {std::max(hi.x, p.x), std::max(hi.y, p.y), std::max(hi.z, p.z)};
    }
    CHECK(norm(lo) < 0.01);
    CHECK(norm(hi - Vec3{2, 1, 1}) < 0.01);
  }
  SECTION("overlapping solids become one hull") {
    const Shell box = make_box({0, 0, 0}, {2, 1, 1}), ball = make_sphere({2, 0.5, 0.5}, 0.8);
    const TriangleMesh m = shrink_wrap(std::vector<Shell>{box, ball}, 0.05);
    check_closed(m);
    CHECK(rel_close(m.volume(), boolean(box, ball, BooleanOp::Union).volume(), 0.01));
  }
  SECTION("a torus keeps its hole") {
    const TriangleMesh m =
        shrink_wrap(std::vector<Shell>{make_torus({0, 0, 0}, {0, 0, 1}, 2.0, 0.5)}, 0.1);
    check_closed(m, 1, 1);
    CHECK(rel_close(m.volume(), 2.0 * pi * pi * 2.0 * 0.25, 0.03));
  }
  SECTION("gaps below the resolution are closed, wider ones are kept") {
    const std::vector<Shell> two{make_box({0, 0, 0}, {1, 1, 1}), make_box({1.04, 0, 0}, {1, 1, 1})};
    const TriangleMesh closed = shrink_wrap(two, 0.05);
    check_closed(closed);
    CHECK(rel_close(closed.volume(), 2.04, 0.01));
    const TriangleMesh apart = shrink_wrap(two, 0.02);
    check_closed(apart, 2);
    CHECK(rel_close(apart.volume(), 2.0, 0.01));
  }
  SECTION("trimmed faces") {
    const Shell part =
        boolean(make_box({-1, -1, -1}, {2, 2, 2}), make_cylinder({0, 0, -2}, {0, 0, 1}, 0.5, 4.0),
                BooleanOp::Difference);
    const TriangleMesh m = shrink_wrap(std::vector<Shell>{part}, 0.05);
    check_closed(m, 1, 1);
    CHECK(rel_close(m.volume(), part.volume(), 0.015));
  }
}

TEST_CASE("shrink_wrap: loose surfaces and errors", "[mesh][shrink_wrap]") {
  const std::vector<TrimmedSurface> sheet{
      TrimmedSurface(make_bilinear({0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0}))};
  ShrinkWrapOptions opt;
  opt.offset = 0.1;
  const TriangleMesh m = shrink_wrap(sheet, 0.05, opt);
  check_closed(m);
  // a plate with rounded rim: 0.2 thick, half cylinders along the sides, sphere octants
  CHECK(rel_close(m.volume(), 0.2 + 4 * 0.5 * pi * 0.01 + 4.0 / 3.0 * pi * 1e-3, 0.03));

  CHECK_THROWS_AS(shrink_wrap(sheet, 0.05), GeometryError);  // no thickness: nothing is left
  CHECK_THROWS_AS(shrink_wrap(sheet, 0.0), GeometryError);
  opt.offset = -0.1;
  CHECK_THROWS_AS(shrink_wrap(sheet, 0.05, opt), GeometryError);
  CHECK_THROWS_AS(shrink_wrap(std::vector<TrimmedSurface>{}, 0.05), GeometryError);
  CHECK_THROWS_AS(shrink_wrap(sheet, 1e-4), GeometryError);  // grid beyond max_nodes
  ShrinkWrapOptions small;
  small.offset = 0.1;
  small.max_nodes = 1000;
  CHECK_THROWS_AS(shrink_wrap(sheet, 0.05, small), GeometryError);
}
