#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <set>
#include <string>
#include <vector>

#include "nurbforge/core/error.hpp"
#include "nurbforge/geometry/transform.hpp"
#include "nurbforge/io/step.hpp"
#include "nurbforge/topology/boolean.hpp"
#include "nurbforge/topology/solids.hpp"

using namespace nf;

namespace {

std::multiset<std::string> labels(const Shell& s) {
  std::multiset<std::string> out;
  for (const auto& f : s.faces()) out.insert(f.label());
  return out;
}

}  // namespace

TEST_CASE("face labels of primitives and simple operations", "[labels]") {
  const auto box = make_box({0, 0, 0}, {2, 2, 1});
  CHECK(labels(box) ==
        std::multiset<std::string>{"bottom", "top", "front", "back", "left", "right"});
  CHECK(labels(make_cylinder({0, 0, 0}, {1, 1, 0}, 1, 2)) ==
        std::multiset<std::string>{"side", "bottom", "top"});
  CHECK(labels(make_cone({0, 0, 0}, {0, 0, 1}, 1, 0, 2)) ==
        std::multiset<std::string>{"side", "bottom"});
  CHECK(labels(make_torus({0, 0, 0}, {0, 0, 1}, 3, 1)) == std::multiset<std::string>{"surface"});
  // flipped, transformed (also mirrored), relabelled
  CHECK(box.faces()[0].flipped().label() == box.faces()[0].label());
  CHECK(labels(box.transformed(Transform::mirror({0, 0, 0}, {1, 0, 0}))) == labels(box));
  std::vector<std::string> names(box.faces().size());
  for (std::size_t i = 0; i < names.size(); ++i) names[i] = "f" + std::to_string(i);
  const auto renamed = box.with_face_labels(names);
  CHECK(renamed.faces_labeled("f3") == std::vector<std::size_t>{3});
  CHECK(renamed.is_closed());
  CHECK_THROWS_AS(box.with_face_labels({"a"}), GeometryError);
}

TEST_CASE("face labels survive Booleans", "[labels]") {
  // plate minus a through cylinder: the plate faces keep their labels, the hole wall is the
  // cylinder's side
  const auto plate = make_box({-2, -2, 0}, {4, 4, 1});
  const auto tool = make_cylinder({0, 0, -1}, {0, 0, 1}, 1, 3);
  const auto cut = boolean(plate, tool, BooleanOp::Difference);
  REQUIRE(cut.is_closed());
  const auto l = labels(cut);
  for (const char* name : {"bottom", "top", "front", "back", "left", "right", "side"})
    CHECK(l.count(name) >= 1);
  CHECK(l.count("") == 0);
  // the hole wall: area 2 pi r h = 2 pi
  double wall = 0;
  for (const auto i : cut.faces_labeled("side")) wall += cut.faces()[i].area();
  CHECK(std::abs(wall - 2 * 3.141592653589793) < 1e-6);

  // union of two labelled boxes: every face of the result comes from a labelled face
  const auto a =
      make_box({0, 0, 0}, {2, 2, 2}).with_face_labels({"a0", "a1", "a2", "a3", "a4", "a5"});
  const auto b =
      make_box({1, 1, 1}, {2, 2, 2}).with_face_labels({"b0", "b1", "b2", "b3", "b4", "b5"});
  const auto u = boolean(a, b, BooleanOp::Union);
  for (const auto& f : u.faces()) CHECK(f.label().size() == 2);
  CHECK(u.faces_labeled("a0").size() == 1);
}

TEST_CASE("face labels round-trip through STEP", "[labels]") {
  const auto c =
      make_cylinder({0, 0, 0}, {0, 0, 1}, 1, 2).with_face_labels({"wall", "inlet", "outlet"});
  const std::string text = to_step({c});
  CHECK(text.find("ADVANCED_FACE('inlet'") != std::string::npos);
  const auto back = from_step(text);
  REQUIRE(back.size() == 1);
  CHECK(labels(back[0]) == std::multiset<std::string>{"wall", "inlet", "outlet"});
}
