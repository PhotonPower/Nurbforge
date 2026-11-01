#include <catch2/catch_test_macros.hpp>

#include <cctype>
#include <cstdio>
#include <fstream>
#include <map>
#include <numbers>
#include <regex>
#include <sstream>
#include <string>
#include <vector>

#include "nurbforge/geometry/primitives.hpp"
#include "nurbforge/io/step.hpp"
#include "nurbforge/modeling/construct.hpp"
#include "nurbforge/topology/boolean.hpp"

using namespace nf;
constexpr double pi = std::numbers::pi;

namespace {

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

// Minimal Part 21 reader: entity instances of the DATA section by id.
struct StepFile {
  std::string header;
  std::map<std::size_t, std::string> entities;  // id -> body (without "#id=" and ";")

  explicit StepFile(const std::string& text) {
    const auto data = text.find("\nDATA;\n");
    REQUIRE(data != std::string::npos);
    header = text.substr(0, data);
    const auto end = text.find("ENDSEC;", data);
    REQUIRE(end != std::string::npos);
    std::string body = text.substr(data + 7, end - data - 7);
    std::string joined;
    for (std::size_t i = 0; i < body.size(); ++i) {
      if (body[i] == '\n' && body.compare(i + 1, 2, "  ") == 0) i += 2;  // continuation line
      else joined += body[i];
    }
    std::istringstream in(joined);
    std::string line;
    while (std::getline(in, line)) {
      REQUIRE(line.size() > 3);
      REQUIRE(line[0] == '#');
      REQUIRE(line.back() == ';');
      const auto eq = line.find('=');
      const std::size_t id = std::stoul(line.substr(1, eq - 1));
      REQUIRE(entities.count(id) == 0);
      entities[id] = line.substr(eq + 1, line.size() - eq - 2);
    }
  }

  // entity type names of an instance (several for complex instances)
  static std::vector<std::string> types(const std::string& body) {
    std::vector<std::string> out;
    if (body[0] != '(') {
      out.push_back(body.substr(0, body.find('(')));
      return out;
    }
    int depth = 0;
    for (std::size_t i = 1; i < body.size(); ++i) {
      if (body[i] == '(') ++depth;
      else if (body[i] == ')') --depth;
      else if (depth == 0 && std::isupper(static_cast<unsigned char>(body[i]))) {
        const auto open = body.find('(', i);
        out.push_back(body.substr(i, open - i));
        i = open - 1;
      }
    }
    return out;
  }

  [[nodiscard]] std::size_t count(const std::string& type) const {
    std::size_t n = 0;
    for (const auto& [id, b] : entities)
      for (const auto& t : types(b))
        if (t == type) ++n;
    return n;
  }

  static std::vector<std::size_t> refs(const std::string& body) {
    std::vector<std::size_t> out;
    static const std::regex r(R"(#(\d+))");
    for (auto it = std::sregex_iterator(body.begin(), body.end(), r); it != std::sregex_iterator(); ++it)
      out.push_back(std::stoul((*it)[1].str()));
    return out;
  }

  // every reference resolves
  [[nodiscard]] bool references_resolve() const {
    for (const auto& [id, b] : entities)
      for (std::size_t r : refs(b))
        if (!entities.count(r)) return false;
    return true;
  }

  // number of ORIENTED_EDGEs referencing each EDGE_CURVE, and the orientation flags
  [[nodiscard]] std::map<std::size_t, std::vector<bool>> edge_uses() const {
    std::map<std::size_t, std::vector<bool>> out;
    for (const auto& [id, b] : entities)
      if (types(b)[0] == "EDGE_CURVE") out[id];
    for (const auto& [id, b] : entities)
      if (types(b)[0] == "ORIENTED_EDGE") out[refs(b).at(0)].push_back(b.find(",.T.)") != std::string::npos);
    return out;
  }
};

void check_common(const StepFile& f) {
  CHECK(f.header.rfind("ISO-10303-21;\nHEADER;\n", 0) == 0);
  CHECK(f.header.find("FILE_SCHEMA(('AUTOMOTIVE_DESIGN") != std::string::npos);
  CHECK(f.references_resolve());
  CHECK(f.count("SHAPE_DEFINITION_REPRESENTATION") == 1);
  CHECK(f.count("PRODUCT") == 1);
  // every real has a decimal point
  static const std::regex bad_real(R"([(,]-?\d+[eE])");
  for (const auto& [id, b] : f.entities) CHECK_FALSE(std::regex_search(b, bad_real));
}

}  // namespace

TEST_CASE("STEP: cube", "[step]") {
  StepOptions opt;
  opt.timestamp = "2026-01-01T00:00:00";
  opt.name = "cube 'one'";
  const std::string text = to_step({box_shell(1.0)}, opt);
  CHECK(text.size() > 1000);
  CHECK(text.substr(text.size() - 18) == "END-ISO-10303-21;\n");
  const StepFile f(text);
  check_common(f);
  CHECK(f.header.find("'cube ''one'''") != std::string::npos);
  CHECK(f.header.find("'2026-01-01T00:00:00'") != std::string::npos);
  CHECK(f.count("MANIFOLD_SOLID_BREP") == 1);
  CHECK(f.count("CLOSED_SHELL") == 1);
  CHECK(f.count("ADVANCED_FACE") == 6);
  CHECK(f.count("FACE_OUTER_BOUND") == 6);
  CHECK(f.count("EDGE_LOOP") == 6);
  CHECK(f.count("EDGE_CURVE") == 12);
  CHECK(f.count("VERTEX_POINT") == 8);
  CHECK(f.count("ORIENTED_EDGE") == 24);
  CHECK(f.count("B_SPLINE_SURFACE_WITH_KNOTS") == 6);
  CHECK(f.count("RATIONAL_B_SPLINE_SURFACE") == 0);
  CHECK(f.count("SI_UNIT") == 3);
  CHECK(text.find("SI_UNIT(.MILLI.,.METRE.)") != std::string::npos);
  for (const auto& [id, uses] : f.edge_uses()) {
    REQUIRE(uses.size() == 2);
    CHECK(uses[0] != uses[1]);  // opposite orientations
  }
  for (const auto& [id, b] : f.entities)
    if (StepFile::types(b)[0] == "ADVANCED_FACE") CHECK(b.substr(b.size() - 5) == ",.T.)");
}

TEST_CASE("STEP: sphere with seam edge and rational surface", "[step]") {
  const StepFile f(to_step({sphere_shell({0, 0, 0}, 2.0)}));
  check_common(f);
  CHECK(f.count("ADVANCED_FACE") == 1);
  CHECK(f.count("RATIONAL_B_SPLINE_SURFACE") == 1);
  CHECK(f.count("BOUNDED_SURFACE") == 1);
  CHECK(f.count("EDGE_CURVE") == 1);    // the seam
  CHECK(f.count("VERTEX_POINT") == 2);  // the poles
  const auto uses = f.edge_uses();
  REQUIRE(uses.size() == 1);
  REQUIRE(uses.begin()->second.size() == 2);
  CHECK(uses.begin()->second[0] != uses.begin()->second[1]);
  for (const auto& [id, b] : f.entities)
    if (StepFile::types(b)[0] == "EDGE_LOOP") CHECK(StepFile::refs(b).size() == 2);  // seam twice
  // complex entity: subtypes in alphabetical order
  bool found = false;
  for (const auto& [id, b] : f.entities)
    if (b[0] == '(' && b.find("RATIONAL_B_SPLINE_SURFACE") != std::string::npos) {
      found = true;
      const std::vector<std::string> expect = {"BOUNDED_SURFACE", "B_SPLINE_SURFACE",
                                               "B_SPLINE_SURFACE_WITH_KNOTS", "GEOMETRIC_REPRESENTATION_ITEM",
                                               "RATIONAL_B_SPLINE_SURFACE", "REPRESENTATION_ITEM", "SURFACE"};
      CHECK(StepFile::types(b) == expect);
    }
  CHECK(found);
}

TEST_CASE("STEP: drilled cube from a Boolean difference", "[step]") {
  const Shell m = boolean(box_shell(1.0), cylinder_shell(0.3, -2.0, 2.0), BooleanOp::Difference);
  REQUIRE(m.is_closed());
  const StepFile f(to_step({m}));
  check_common(f);
  CHECK(f.count("MANIFOLD_SOLID_BREP") == 1);
  CHECK(f.count("ADVANCED_FACE") == m.faces().size());
  CHECK(f.count("EDGE_CURVE") == m.edges().size());
  CHECK(f.count("FACE_BOUND") == 2);  // the holes in the top and bottom faces
  for (const auto& [id, uses] : f.edge_uses()) {
    REQUIRE(uses.size() == 2);
    CHECK(uses[0] != uses[1]);
  }
}

TEST_CASE("STEP: void, open shell, units, file output", "[step]") {
  SECTION("cube with a spherical void") {
    const Shell m = boolean(box_shell(1.0), sphere_shell({0.1, -0.1, 0.2}, 0.5), BooleanOp::Difference);
    const StepFile f(to_step({m}));
    check_common(f);
    CHECK(f.count("BREP_WITH_VOIDS") == 1);
    CHECK(f.count("ORIENTED_CLOSED_SHELL") == 1);
    CHECK(f.count("CLOSED_SHELL") == 2);
    CHECK(f.count("MANIFOLD_SOLID_BREP") == 0);
  }
  SECTION("open shell and a solid") {
    const Shell plate = Shell::sew({TrimmedSurface(make_rectangle({0, 0, 0}, {1, 0, 0}, {0, 1, 0}))});
    StepOptions opt;
    opt.unit = StepLengthUnit::Inch;
    const std::string text = to_step({plate, box_shell(1.0)}, opt);
    const StepFile f(text);
    check_common(f);
    CHECK(f.count("OPEN_SHELL") == 1);
    CHECK(f.count("SHELL_BASED_SURFACE_MODEL") == 1);
    CHECK(f.count("SHAPE_REPRESENTATION") == 1);  // both kinds of items in one representation
    CHECK(f.count("MANIFOLD_SOLID_BREP") == 1);
    CHECK(f.count("CONVERSION_BASED_UNIT") == 1);
    const StepFile open_only(to_step({plate}));
    CHECK(open_only.count("MANIFOLD_SURFACE_SHAPE_REPRESENTATION") == 1);
    CHECK(open_only.count("ADVANCED_BREP_SHAPE_REPRESENTATION") == 0);
  }
  SECTION("inside-out sphere is written with reversed face sense") {
    const Shell inv = Shell::sew({TrimmedSurface(
        revolve(make_arc({0, 0, 0}, {0, 0, -1}, {1, 0, 0}, 1.0, 0.0, pi), {0, 0, 0}, {0, 0, -1}, 2 * pi))});
    REQUIRE(inv.volume(1e-6) < 0.0);
    const StepFile f(to_step({inv}));
    for (const auto& [id, b] : f.entities) {
      const auto t = StepFile::types(b)[0];
      if (t == "ADVANCED_FACE" || t == "FACE_OUTER_BOUND") CHECK(b.substr(b.size() - 5) == ",.F.)");
    }
  }
  SECTION("write_step") {
    const std::string path = "nf_test_step_output.step";
    write_step(path, {box_shell(0.5)});
    std::ifstream in(path);
    std::stringstream ss;
    ss << in.rdbuf();
    in.close();
    std::remove(path.c_str());
    CHECK(StepFile(ss.str()).count("ADVANCED_FACE") == 6);
    CHECK_THROWS(write_step("no_such_dir/x/y.step", {box_shell(0.5)}));
  }
}
