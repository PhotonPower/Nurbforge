#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <numbers>
#include <vector>

#include "nurbforge/core/error.hpp"
#include "nurbforge/geometry/primitives.hpp"
#include "nurbforge/geometry/transform.hpp"
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

// ------------------------------------------------------------------ bodies with several shells

namespace {
std::vector<TrimmedSurface> faces_of(const Shell& s, bool flip) {
  std::vector<TrimmedSurface> out;
  for (const auto& f : s.faces()) out.push_back(flip ? f.flipped() : f);
  return out;
}
Shell hollow_cube(double void_r, bool flip_void) {  // [-1, 1]^3 with a spherical void
  auto faces = faces_of(box_shell(1.0), false);
  for (auto& f : faces_of(sphere_shell({0, 0, 0}, void_r), flip_void)) faces.push_back(std::move(f));
  return Shell::sew(std::move(faces));
}
}  // namespace

TEST_CASE("orientation of shells with several components (voids, lumps)", "[shell][bodies]") {
  const double r = 0.5, vs = 4.0 / 3.0 * pi * r * r * r;
  for (bool flip_void : {false, true}) {
    const Shell h = hollow_cube(r, flip_void).oriented();
    REQUIRE(h.is_closed());
    CHECK(rel_close(h.volume(1e-8), 8.0 - vs, 1e-7));
    const auto comps = h.components();
    REQUIRE(comps.size() == 2);
    const double v0 = comps[0].volume(1e-8), v1 = comps[1].volume(1e-8);
    CHECK(rel_close(std::min(v0, v1), -vs, 1e-7));  // the void bounds negative volume
    CHECK(rel_close(std::max(v0, v1), 8.0, 1e-7));
  }
  // two separate spheres, one given inside out: both become outward
  auto faces = faces_of(sphere_shell({0, 0, 0}, 1.0), false);
  for (auto& f : faces_of(sphere_shell({3, 0, 0}, 1.0), true)) faces.push_back(std::move(f));
  const Shell two = Shell::sew(std::move(faces)).oriented();
  CHECK(two.components().size() == 2);
  CHECK(rel_close(two.volume(1e-8), 2.0 * 4.0 / 3.0 * pi, 1e-7));
}

TEST_CASE("Booleans with hollow inputs and several lumps", "[boolean][bodies]") {
  const double r = 0.5, vs = 4.0 / 3.0 * pi * r * r * r;
  const Shell hollow = hollow_cube(r, true);
  SECTION("union with a ball floating in the void") {
    const double q = 0.25;
    const Shell u = boolean(hollow, sphere_shell({0, 0, 0}, q), BooleanOp::Union);
    REQUIRE(u.is_closed());
    CHECK(u.components().size() == 3);
    CHECK(rel_close(u.volume(1e-8), 8.0 - vs + 4.0 / 3.0 * pi * q * q * q, 1e-6));
  }
  SECTION("intersection cutting through the void") {
    const Shell i = boolean(hollow, box_shell({0, -2, -2}, {2, 2, 2}), BooleanOp::Intersection);
    REQUIRE(i.is_closed());
    CHECK(rel_close(i.volume(1e-8), 4.0 - 0.5 * vs, 1e-6));
  }
  SECTION("difference splitting a cube into two lumps") {
    const Shell d = boolean(box_shell(1.0), box_shell({-0.25, -2, -2}, {0.25, 2, 2}), BooleanOp::Difference);
    REQUIRE(d.is_closed());
    CHECK(d.components().size() == 2);
    CHECK(rel_close(d.volume(1e-8), 6.0, 1e-6));
    CHECK(rel_close(d.oriented().volume(1e-8), 6.0, 1e-6));
  }
}

TEST_CASE("Boolean with a plane through both poles of a sphere", "[boolean][bodies]") {
  // the intersection curves are meridians ending at the poles (degenerate boundary): the marching
  // must neither reverse there nor create knot-line points between unrelated pole parameters
  const double r = 0.5;
  const Shell half = boolean(sphere_shell({0, 0, 0}, r), box_shell({0, -2, -2}, {2, 2, 2}), BooleanOp::Intersection);
  REQUIRE(half.is_closed());
  CHECK(half.components().size() == 1);
  CHECK(rel_close(half.volume(1e-8), 2.0 / 3.0 * pi * r * r * r, 1e-7));
}

// ------------------------------------------------------------------ robustness (stress-test finds)

namespace {
Shell offset_cylinder(const Vec3& base, double r, double h) {  // axis z
  const auto side = extrude(make_circle(base, {0, 0, 1}, r), {0, 0, h});
  const auto lid = [&](double z) {
    return TrimmedSurface(make_rectangle({base.x - r, base.y - r, z}, {2 * r, 0, 0}, {0, 2 * r, 0}),
                          circle_loop(0.5, 0.5, 0.5));
  };
  return Shell::sew({TrimmedSurface(side), lid(base.z), lid(base.z + h)});
}
// closed, consistently oriented results and V(A u B) + V(A n B) = V(A) + V(B), V(A - B) = V(A) - V(A n B)
void check_identities(const Shell& a, const Shell& b) {
  const double va = a.oriented().volume(1e-7), vb = b.oriented().volume(1e-7);
  const Shell u = boolean(a, b, BooleanOp::Union), i = boolean(a, b, BooleanOp::Intersection),
              d = boolean(a, b, BooleanOp::Difference);
  for (const Shell* r : {&u, &i, &d}) {
    REQUIRE(r->is_closed());
    CHECK(r->is_consistently_oriented());
  }
  const double vu = u.volume(1e-7), vi = i.volume(1e-7), vd = d.volume(1e-7);
  CHECK(vi > 0.0);
  CHECK(rel_close(vu + vi, va + vb, 1e-5));
  CHECK(rel_close(vd, va - vi, 1e-5));
}
}  // namespace

TEST_CASE("Booleans: intersection curves ending on a trim circle", "[boolean][robust]") {
  // box faces cut a cylinder lid: the cut ends lie on the lid's trim circle, which the polygon
  // test (64-gon) misplaced by up to 6e-4 – the cuts ended short and the lid stayed whole
  check_identities(cylinder_shell(0.7, 0.0, 2.0),
                   box_shell({-0.97734720018362475, -0.67524846625610246, -0.051429789902464274},
                             {0.33479161330348473, 0.39412143103613317, 0.72290733704003352}));
}

TEST_CASE("Booleans: closest point near a pole", "[boolean][robust]") {
  // projecting cylinder points onto the sphere near its pole split degenerate patches below the
  // parametric resolution (closest_point threw)
  check_identities(sphere_shell({0, 0, 0}, 1.199961155503575),
                   offset_cylinder({-0.40259733731404257, -0.46281045397980786, -1.8177877666917477},
                                   0.35243433601614005, 1.380810480634058));
}

TEST_CASE("Booleans: long marching steps on small faces", "[boolean][robust]") {
  // a straight intersection line on a small box face: a step covering more than a quarter of the
  // face's parameter range was taken for a jump over a pole and the curve was cut short
  check_identities(cylinder_shell(0.7, 0.0, 2.0),
                   box_shell({-0.73456574283628506, -0.27082589505852428, 0.21489243439744388},
                             {-0.069787101918596051, 0.2784369881762836, 0.67397747262179897}));
}

// ------------------------------------------------------------------ tangential contact

namespace {
// expected volumes of union, intersection (0: empty result) and difference; every non-empty
// result closed and consistently oriented
void check_contact(const Shell& a, const Shell& b, double vu, double vi, double vd) {
  const double want[3] = {vu, vi, vd};
  const BooleanOp ops[3] = {BooleanOp::Union, BooleanOp::Intersection, BooleanOp::Difference};
  for (int k = 0; k < 3; ++k) {
    const Shell r = boolean(a, b, ops[k]);
    if (want[k] == 0.0) {
      CHECK(r.faces().empty());
      continue;
    }
    REQUIRE(r.is_closed());
    CHECK(r.is_consistently_oriented());
    CHECK(rel_close(r.volume(1e-7), want[k], 1e-6));
  }
}
constexpr double unit_ball = 4.0 / 3.0 * pi;
}  // namespace

TEST_CASE("Booleans: spheres touching at a point", "[boolean][tangent]") {
  // externally: two lumps / nothing / A
  check_contact(sphere_shell({0, 0, 0}, 1.0), sphere_shell({1.5, 0, 0}, 0.5), unit_ball * 1.125, 0.0, unit_ball);
  // internally, contact on the seam, off the seam and at a pole (the marcher produced
  // micrometre-long curves around a contact on the seam)
  const double r3 = 0.064;  // 0.4^3
  for (const Vec3& c : {Vec3{0.6, 0, 0}, Vec3{0.6 * std::cos(1.0), 0.6 * std::sin(1.0), 0}, Vec3{0, 0, 0.6}})
    check_contact(sphere_shell({0, 0, 0}, 1.0), sphere_shell(c, 0.4), unit_ball, unit_ball * r3, unit_ball * (1 - r3));
}

TEST_CASE("Booleans: sphere touching a box face", "[boolean][tangent]") {
  check_contact(sphere_shell({0, 0, 0}, 1.0), box_shell({-1, -1, 1}, {1, 1, 2}), unit_ball + 4.0, 0.0, unit_ball);
  check_contact(box_shell({-1, -1, -1}, {1, 1, 1.5}), sphere_shell({0, 0, 0}, 1.0), 10.0, unit_ball, 10.0 - unit_ball);
}

TEST_CASE("Booleans: cylinder touching a box face along a line", "[boolean][tangent]") {
  // the face's widest interior point lies on the contact line, where ray casting is undecided
  const double cyl = pi * 0.49 * 2.0;
  check_contact(cylinder_shell(0.7, 0.0, 2.0), box_shell({0.7, -1, 0.5}, {1.5, 1, 1.5}), cyl + 1.6, 0.0, cyl);
  check_contact(cylinder_shell(0.7, 0.0, 2.0), box_shell({-1, 0.7, 0.5}, {1, 1.5, 1.5}), cyl + 1.6, 0.0, cyl);
}

TEST_CASE("Booleans: cylinder inside a box touching two sides", "[boolean][tangent]") {
  // the cut circles on the box faces touch the face boundaries: the intersection curve must not
  // end there, and the band of the cylinder between them lies on the contact lines
  const double box = 1.7 * 1.7 * 1.0, mid = pi * 0.49;
  check_contact(box_shell({-0.7, -0.7, 0.5}, {1, 1, 1.5}), cylinder_shell(0.7, 0.0, 2.0),
                box + 2.0 * pi * 0.49 - mid, mid, box - mid);
}

TEST_CASE("Booleans: cylinders touching along a line", "[boolean][tangent]") {
  // curved surfaces touching along a line: the marcher produced scattered fragments
  const double a = pi * 0.49 * 2.0, b = pi * 0.25 * 2.0;
  check_contact(cylinder_shell(0.7, 0.0, 2.0),
                offset_cylinder({1.2 * std::cos(0.4), 1.2 * std::sin(0.4), 0.5}, 0.5, 2.0), a + b, 0.0, a);
}

namespace {
Shell rotated_z(const Shell& s, double angle) {
  std::vector<TrimmedSurface> faces;
  for (const auto& f : s.faces())
    faces.emplace_back(f.surface().transformed(Transform::rotation({0, 0, 0}, {0, 0, 1}, angle)), f.outer(), f.holes());
  return Shell::sew(std::move(faces));
}
}  // namespace

TEST_CASE("Booleans: tilted box face touching a cylinder", "[boolean][tangent]") {
  // the patch boundaries of the cylinder next to the contact line run just off the tilted plane:
  // the curve/surface intersection subdivided them down to the size limit (minutes)
  const double r = 0.7021256186899939, z0 = -0.21300335723297636, h = 1.2982978724505516, w = 0.67037666446677946;
  const Shell b = rotated_z(box_shell({r, -1.0, z0}, {r + w, 1.0, z0 + h}), 5.1897199994457743);
  const double cyl = pi * r * r * 2.0;
  check_contact(cylinder_shell(r, 0.0, 2.0), b, cyl + 2.0 * w * h, 0.0, cyl);
}

TEST_CASE("Booleans: cylinder inside a tilted box touching two sides", "[boolean][tangent]") {
  // the two ends of the intersection curve at a touching point overlap slightly (the position
  // along the tangent is only determined to ~sqrt(2 r tol)) and must still be joined
  const double r = 0.62334641101607957, z0 = 0.66660399501980083, ex = 0.24319203551411819, ey = 0.47548877062290407;
  const Shell a = rotated_z(box_shell({-r, -r, z0}, {r + ex, r + ey, z0 + 0.8}), 4.8653577882165431);
  const double box = (2 * r + ex) * (2 * r + ey) * 0.8, mid = pi * r * r * 0.8;
  check_contact(a, cylinder_shell(r, 0.0, 2.0), box + pi * r * r * 2.0 - mid, mid, box - mid);
}
