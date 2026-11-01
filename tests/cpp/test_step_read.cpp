#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <numbers>
#include <string>
#include <vector>

#include "nurbforge/core/error.hpp"
#include "nurbforge/geometry/primitives.hpp"
#include "nurbforge/io/step.hpp"
#include "nurbforge/modeling/construct.hpp"
#include "nurbforge/topology/boolean.hpp"

using namespace nf;
using Catch::Approx;
constexpr double pi = std::numbers::pi;

namespace {

Shell sphere_shell(const Vec3& c, double R) {
  return Shell::sew({TrimmedSurface(revolve(make_arc(c, {0, 0, -1}, {1, 0, 0}, R, 0.0, pi), c, {0, 0, 1}, 2 * pi))});
}

std::vector<TrimmedSurface> box_faces(double h) {  // cube [-h, h]^3
  const double s = 2 * h;
  return {TrimmedSurface(make_rectangle({-h, -h, -h}, {s, 0, 0}, {0, s, 0})),
          TrimmedSurface(make_rectangle({-h, -h, -h}, {s, 0, 0}, {0, 0, s})),
          TrimmedSurface(make_rectangle({-h, h, -h}, {s, 0, 0}, {0, 0, s})),
          TrimmedSurface(make_rectangle({-h, -h, -h}, {0, s, 0}, {0, 0, s})),
          TrimmedSurface(make_rectangle({h, -h, -h}, {0, s, 0}, {0, 0, s})),
          TrimmedSurface(make_rectangle({-h, -h, h}, {s, 0, 0}, {0, s, 0}))};
}

Shell box_shell(double h) { return Shell::sew(box_faces(h)); }

Shell cylinder_shell(double R, double z0, double z1) {  // axis z
  const auto side = extrude(make_circle({0, 0, z0}, {0, 0, 1}, R), {0, 0, z1 - z0});
  const auto lid = [&](double z) {
    return TrimmedSurface(make_rectangle({-R, -R, z}, {2 * R, 0, 0}, {0, 2 * R, 0}), circle_loop(0.5, 0.5, 0.5));
  };
  return Shell::sew({TrimmedSurface(side), lid(z0), lid(z1)});
}

constexpr double kRel = 1e-7;  // integration accuracy of the property comparisons

double abs_volume(const Shell& s) {
  return std::abs(s.is_consistently_oriented() ? s.volume(kRel) : s.oriented().volume(kRel));
}

/// Writes and reads back; checks one shell with the same face count, closedness, volume and area.
Shell roundtrip(const Shell& original) {
  const auto shells = from_step(to_step({original}));
  REQUIRE(shells.size() == 1);
  const Shell& s = shells[0];
  CHECK(s.faces().size() == original.faces().size());
  CHECK(s.is_closed() == original.is_closed());
  CHECK(s.is_consistently_oriented());
  CHECK(s.area(kRel) == Approx(original.area(kRel)).epsilon(1e-6));
  if (original.is_closed()) {
    const double v = s.volume(kRel);
    CHECK(v > 0.0);
    CHECK(v == Approx(abs_volume(original)).epsilon(1e-6));
  }
  return s;
}

std::string wrap(const std::string& data, const std::string& unit = "(LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT(.MILLI.,.METRE.))") {
  return "ISO-10303-21;\nHEADER;\nFILE_DESCRIPTION(('test'),'2;1');\nFILE_NAME('t','2026-01-01T00:00:00',(''),(''),'',"
         "'','');\nFILE_SCHEMA(('AUTOMOTIVE_DESIGN { 1 0 10303 214 1 1 1 1 }'));\nENDSEC;\nDATA;\n" +
         data + "#900=" + unit +
         ";\n#901=(GEOMETRIC_REPRESENTATION_CONTEXT(3) GLOBAL_UNIT_ASSIGNED_CONTEXT((#900)) "
         "REPRESENTATION_CONTEXT('',''));\nENDSEC;\nEND-ISO-10303-21;\n";
}

// Unit cube [0, 1]^3 from PLANE / LINE entities (with comments, odd line breaks, a FACE_BOUND
// with orientation .F. and a face with same_sense .F.).
const std::string kCube = R"(
#1=CARTESIAN_POINT('',(0.,0.,0.));#2=CARTESIAN_POINT('',(1.,0.,0.));
#3=CARTESIAN_POINT('',(1.,1.,0.));#4=CARTESIAN_POINT('',(0.,1.,0.));
#5=CARTESIAN_POINT('',(0.,0.,1.));#6=CARTESIAN_POINT('',(1.,0.,1.));
#7=CARTESIAN_POINT('',(1.,1.,1.));#8=CARTESIAN_POINT('',(0.,1.,1.));
/* vertices */
#11=VERTEX_POINT('',#1);#12=VERTEX_POINT('',#2);#13=VERTEX_POINT('',#3);#14=VERTEX_POINT('',#4);
#15=VERTEX_POINT('',#5);#16=VERTEX_POINT('',#6);#17=VERTEX_POINT('',#7);#18=VERTEX_POINT('',#8);
#20=DIRECTION('',(1.,0.,0.));#21=DIRECTION('',(0.,1.,0.));#22=DIRECTION('',(0.,0.,1.));
#23=DIRECTION('',(-1.,0.,0.));#24=DIRECTION('',(0.,-1.,0.));
#25=VECTOR('',#20,1.);#26=VECTOR('',#21,1.);#27=VECTOR('',#22,
  1.);
#31=LINE('',#1,#25);#32=LINE('',#2,#26);#33=LINE('',#4,#25);#34=LINE('',#1,#26);
#35=LINE('',#5,#25);#36=LINE('',#6,#26);#37=LINE('',#8,#25);#38=LINE('',#5,#26);
#39=LINE('',#1,#27);#40=LINE('',#2,#27);#41=LINE('',#3,#27);#42=LINE('',#4,#27);
#51=EDGE_CURVE('',#11,#12,#31,.T.);#52=EDGE_CURVE('',#12,#13,#32,.T.);
#53=EDGE_CURVE('',#14,#13,#33,.T.);#54=EDGE_CURVE('',#11,#14,#34,.T.);
#55=EDGE_CURVE('',#15,#16,#35,.T.);#56=EDGE_CURVE('',#16,#17,#36,.T.);
#57=EDGE_CURVE('',#18,#17,#37,.T.);#58=EDGE_CURVE('',#15,#18,#38,.T.);
#59=EDGE_CURVE('',#11,#15,#39,.T.);#60=EDGE_CURVE('',#12,#16,#40,.T.);
#61=EDGE_CURVE('',#13,#17,#41,.T.);#62=EDGE_CURVE('',#14,#18,#42,.T.);
#100=AXIS2_PLACEMENT_3D('',#1,#22,#20);#101=PLANE('',#100);
#102=ORIENTED_EDGE('',*,*,#54,.T.);#103=ORIENTED_EDGE('',*,*,#53,.T.);
#104=ORIENTED_EDGE('',*,*,#52,.F.);#105=ORIENTED_EDGE('',*,*,#51,.F.);
#106=EDGE_LOOP('',(#102,#103,#104,#105));#107=FACE_OUTER_BOUND('',#106,.T.);
#108=ADVANCED_FACE('bottom',(#107),#101,.F.);
#110=AXIS2_PLACEMENT_3D('',#5,#22,#20);#111=PLANE('',#110);
#112=ORIENTED_EDGE('',*,*,#55,.T.);#113=ORIENTED_EDGE('',*,*,#56,.T.);
#114=ORIENTED_EDGE('',*,*,#57,.F.);#115=ORIENTED_EDGE('',*,*,#58,.F.);
#116=EDGE_LOOP('',(#112,#113,#114,#115));#117=FACE_OUTER_BOUND('',#116,.T.);
#118=ADVANCED_FACE('top',(#117),#111,.T.);
#120=AXIS2_PLACEMENT_3D('',#1,#24,#20);#121=PLANE('',#120);
#122=ORIENTED_EDGE('',*,*,#51,.T.);#123=ORIENTED_EDGE('',*,*,#60,.T.);
#124=ORIENTED_EDGE('',*,*,#55,.F.);#125=ORIENTED_EDGE('',*,*,#59,.F.);
#126=EDGE_LOOP('',(#122,#123,#124,#125));#127=FACE_BOUND('',#126,.T.);
#128=ADVANCED_FACE('front',(#127),#121,.T.);
#130=AXIS2_PLACEMENT_3D('',#4,#21,#20);#131=PLANE('',#130);
#132=ORIENTED_EDGE('',*,*,#62,.T.);#133=ORIENTED_EDGE('',*,*,#57,.T.);
#134=ORIENTED_EDGE('',*,*,#61,.F.);#135=ORIENTED_EDGE('',*,*,#53,.F.);
#136=EDGE_LOOP('',(#132,#133,#134,#135));#137=FACE_OUTER_BOUND('',#136,.T.);
#138=ADVANCED_FACE('back',(#137),#131,.T.);
#140=AXIS2_PLACEMENT_3D('',#1,#23,#21);#141=PLANE('',#140);
#142=ORIENTED_EDGE('',*,*,#59,.T.);#143=ORIENTED_EDGE('',*,*,#58,.T.);
#144=ORIENTED_EDGE('',*,*,#62,.F.);#145=ORIENTED_EDGE('',*,*,#54,.F.);
#146=EDGE_LOOP('',(#142,#143,#144,#145));#147=FACE_OUTER_BOUND('',#146,.T.);
#148=ADVANCED_FACE('left',(#147),#141,.T.);
#150=AXIS2_PLACEMENT_3D('',#2,#20,#21);#151=PLANE('',#150);
#152=ORIENTED_EDGE('',*,*,#60,.T.);#153=ORIENTED_EDGE('',*,*,#56,.T.);
#154=ORIENTED_EDGE('',*,*,#61,.F.);#155=ORIENTED_EDGE('',*,*,#52,.F.);
#156=EDGE_LOOP('',(#152,#153,#154,#155));#157=FACE_BOUND('',#156,.F.); /* reversed loop */
#158=ADVANCED_FACE('right',(#157),#151,.T.);
#160=CLOSED_SHELL('',(#108,#118,#128,#138,#148,#158));
#161=MANIFOLD_SOLID_BREP('it''s a cube',#160);
)";

// Cylinder R = 0.5, z in [0, 2] from CYLINDRICAL_SURFACE, CIRCLE and a seam LINE.
const std::string kCylinder = R"(
#1=CARTESIAN_POINT('',(0.,0.,0.));#2=CARTESIAN_POINT('',(0.,0.,2.));
#3=CARTESIAN_POINT('',(0.5,0.,0.));#4=CARTESIAN_POINT('',(0.5,0.,2.));
#5=DIRECTION('',(0.,0.,1.));#6=DIRECTION('',(1.,0.,0.));
#7=AXIS2_PLACEMENT_3D('',#1,#5,#6);#8=AXIS2_PLACEMENT_3D('',#2,#5,#6);
#10=VERTEX_POINT('',#3);#11=VERTEX_POINT('',#4);
#12=CIRCLE('',#7,0.5);#13=CIRCLE('',#8,0.5);#14=VECTOR('',#5,1.);#15=LINE('',#3,#14);
#20=EDGE_CURVE('',#10,#10,#12,.T.);#21=EDGE_CURVE('',#11,#11,#13,.T.);
#22=EDGE_CURVE('',#10,#11,#15,.T.);
#30=CYLINDRICAL_SURFACE('',#7,0.5);
#31=ORIENTED_EDGE('',*,*,#20,.T.);#32=ORIENTED_EDGE('',*,*,#22,.T.);
#33=ORIENTED_EDGE('',*,*,#21,.F.);#34=ORIENTED_EDGE('',*,*,#22,.F.);
#35=EDGE_LOOP('',(#31,#32,#33,#34));#36=FACE_BOUND('',#35,.T.);
#37=ADVANCED_FACE('',(#36),#30,.T.);
#40=PLANE('',#7);#41=ORIENTED_EDGE('',*,*,#20,.F.);#42=EDGE_LOOP('',(#41));
#43=FACE_BOUND('',#42,.T.);#44=ADVANCED_FACE('',(#43),#40,.F.);
#50=PLANE('',#8);#51=ORIENTED_EDGE('',*,*,#21,.T.);#52=EDGE_LOOP('',(#51));
#53=FACE_BOUND('',#52,.T.);#54=ADVANCED_FACE('',(#53),#50,.T.);
#60=CLOSED_SHELL('',(#37,#44,#54));#61=MANIFOLD_SOLID_BREP('',#60);
)";

// Half cylinder R = 1, y >= 0, z in [0, 1]: partial cylindrical face bounded by two arcs.
const std::string kHalfCylinder = R"(
#1=CARTESIAN_POINT('',(0.,0.,0.));#2=CARTESIAN_POINT('',(0.,0.,1.));
#3=CARTESIAN_POINT('',(1.,0.,0.));#4=CARTESIAN_POINT('',(-1.,0.,0.));
#5=CARTESIAN_POINT('',(1.,0.,1.));#6=CARTESIAN_POINT('',(-1.,0.,1.));
#7=DIRECTION('',(0.,0.,1.));#8=DIRECTION('',(1.,0.,0.));#9=DIRECTION('',(-1.,0.,0.));
#10=DIRECTION('',(0.,-1.,0.));
#11=AXIS2_PLACEMENT_3D('',#1,#7,#8);#12=AXIS2_PLACEMENT_3D('',#2,#7,#8);
#13=AXIS2_PLACEMENT_3D('',#1,#10,#8);
#20=VERTEX_POINT('',#3);#21=VERTEX_POINT('',#4);#22=VERTEX_POINT('',#5);#23=VERTEX_POINT('',#6);
#24=VECTOR('',#7,1.);#25=VECTOR('',#9,2.);
#30=EDGE_CURVE('',#20,#21,#31,.T.);#31=CIRCLE('',#11,1.);
#32=EDGE_CURVE('',#22,#23,#33,.T.);#33=CIRCLE('',#12,1.);
#34=EDGE_CURVE('',#20,#22,#35,.T.);#35=LINE('',#3,#24);
#36=EDGE_CURVE('',#21,#23,#37,.T.);#37=LINE('',#4,#24);
#38=EDGE_CURVE('',#20,#21,#39,.T.);#39=LINE('',#3,#25);
#40=EDGE_CURVE('',#22,#23,#41,.T.);#41=LINE('',#5,#25);
#50=CYLINDRICAL_SURFACE('',#11,1.);
#51=EDGE_LOOP('',(#52,#53,#54,#55));
#52=ORIENTED_EDGE('',*,*,#30,.T.);#53=ORIENTED_EDGE('',*,*,#36,.T.);
#54=ORIENTED_EDGE('',*,*,#32,.F.);#55=ORIENTED_EDGE('',*,*,#34,.F.);
#56=ADVANCED_FACE('',(#57),#50,.T.);#57=FACE_OUTER_BOUND('',#51,.T.);
#60=PLANE('',#13);
#61=EDGE_LOOP('',(#62,#63,#64,#65));
#62=ORIENTED_EDGE('',*,*,#34,.T.);#63=ORIENTED_EDGE('',*,*,#40,.T.);
#64=ORIENTED_EDGE('',*,*,#36,.F.);#65=ORIENTED_EDGE('',*,*,#38,.F.);
#66=ADVANCED_FACE('',(#67),#60,.T.);#67=FACE_OUTER_BOUND('',#61,.T.);
#70=PLANE('',#11);
#71=EDGE_LOOP('',(#72,#73));#72=ORIENTED_EDGE('',*,*,#30,.F.);#73=ORIENTED_EDGE('',*,*,#38,.T.);
#74=ADVANCED_FACE('',(#75),#70,.F.);#75=FACE_OUTER_BOUND('',#71,.T.);
#80=PLANE('',#12);
#81=EDGE_LOOP('',(#82,#83));#82=ORIENTED_EDGE('',*,*,#32,.T.);#83=ORIENTED_EDGE('',*,*,#40,.F.);
#84=ADVANCED_FACE('',(#85),#80,.T.);#85=FACE_OUTER_BOUND('',#81,.T.);
#90=CLOSED_SHELL('',(#56,#66,#74,#84));#91=MANIFOLD_SOLID_BREP('',#90);
)";

// Sphere R = 2 about (1, 2, 3) as OpenCASCADE writes it: one face bounded by a VERTEX_LOOP.
const std::string kSphere = R"(
#1=CARTESIAN_POINT('',(1.,2.,3.));#2=DIRECTION('',(0.,0.,1.));#3=DIRECTION('',(1.,0.,0.));
#4=AXIS2_PLACEMENT_3D('',#1,#2,#3);#5=SPHERICAL_SURFACE('',#4,2.);
#6=CARTESIAN_POINT('',(1.,2.,1.));#7=VERTEX_POINT('',#6);#8=VERTEX_LOOP('',#7);
#9=FACE_BOUND('',#8,.T.);#10=ADVANCED_FACE('',(#9),#5,.T.);
#11=CLOSED_SHELL('',(#10));#12=MANIFOLD_SOLID_BREP('',#11);
)";

// Hemisphere R = 1, z >= 0: spherical cap bounded by the equator, and a disc.
const std::string kHemisphere = R"(
#1=CARTESIAN_POINT('',(0.,0.,0.));#2=DIRECTION('',(0.,0.,1.));#3=DIRECTION('',(1.,0.,0.));
#4=AXIS2_PLACEMENT_3D('',#1,#2,#3);#5=SPHERICAL_SURFACE('',#4,1.);#6=PLANE('',#4);
#7=CARTESIAN_POINT('',(0.,-1.,0.));#8=VERTEX_POINT('',#7);#9=CIRCLE('',#4,1.);
#10=EDGE_CURVE('',#8,#8,#9,.T.);
#11=ORIENTED_EDGE('',*,*,#10,.T.);#12=EDGE_LOOP('',(#11));#13=FACE_BOUND('',#12,.T.);
#14=ADVANCED_FACE('',(#13),#5,.T.);
#15=ORIENTED_EDGE('',*,*,#10,.F.);#16=EDGE_LOOP('',(#15));#17=FACE_BOUND('',#16,.T.);
#18=ADVANCED_FACE('',(#17),#6,.F.);
#19=CLOSED_SHELL('',(#14,#18));#20=MANIFOLD_SOLID_BREP('',#19);
)";

Shell read_one(const std::string& text, const StepReadOptions& opt = {}) {
  const auto shells = from_step(text, opt);
  REQUIRE(shells.size() == 1);
  return shells[0];
}

}  // namespace

TEST_CASE("STEP read: round trip of a cube", "[step_read]") {
  const Shell s = roundtrip(box_shell(1.0));
  CHECK(s.volume(kRel) == Approx(8.0).epsilon(1e-9));
  CHECK(s.edges().size() == 12);
}

TEST_CASE("STEP read: round trip of a sphere", "[step_read]") {
  const Shell s = roundtrip(sphere_shell({0.5, -1, 2}, 1.5));
  CHECK(s.volume(kRel) == Approx(4.0 / 3.0 * pi * 1.5 * 1.5 * 1.5).epsilon(1e-7));
}

TEST_CASE("STEP read: round trip of a cylinder with lids", "[step_read]") {
  const Shell s = roundtrip(cylinder_shell(0.7, -1.0, 2.0));
  CHECK(s.volume(kRel) == Approx(pi * 0.49 * 3.0).epsilon(1e-7));
}

TEST_CASE("STEP read: round trip of a Boolean union", "[step_read]") {
  const Shell u = boolean(sphere_shell({0, 0, 0}, 1.0), sphere_shell({1, 0, 0}, 1.0), BooleanOp::Union);
  roundtrip(u);
}

TEST_CASE("STEP read: round trip of a solid with a void", "[step_read]") {
  const Shell v = boolean(box_shell(1.0), sphere_shell({0.1, -0.1, 0.2}, 0.5), BooleanOp::Difference);
  const std::string text = to_step({v});
  CHECK(text.find("BREP_WITH_VOIDS") != std::string::npos);
  const Shell s = roundtrip(v);
  CHECK(s.components().size() == 2);
  CHECK(s.volume(kRel) == Approx(8.0 - 4.0 / 3.0 * pi * 0.125).epsilon(1e-7));
}

TEST_CASE("STEP read: round trip of an open shell", "[step_read]") {
  auto faces = box_faces(1.0);
  faces.pop_back();  // no lid
  faces[0] = faces[0].with_hole(circle_loop(0.5, 0.5, 0.25));
  const Shell open = Shell::sew(std::move(faces));
  REQUIRE_FALSE(open.is_closed());
  const Shell s = roundtrip(open);
  CHECK(s.free_edges().size() == open.free_edges().size());
  CHECK(s.area() == Approx(20.0 - pi * 0.25).epsilon(1e-7));
}

TEST_CASE("STEP read: several shells and the writer's options", "[step_read]") {
  StepOptions w;
  w.unit = StepLengthUnit::Metre;
  const auto shells = from_step(to_step({box_shell(1.0), sphere_shell({5, 0, 0}, 1.0)}, w));
  REQUIRE(shells.size() == 2);
  CHECK(shells[0].volume(kRel) == Approx(8e9).epsilon(1e-9));  // m^3 -> mm^3
  CHECK(shells[1].volume(kRel) == Approx(4.0 / 3.0 * pi * 1e9).epsilon(1e-7));
  StepReadOptions r;
  r.unit = StepLengthUnit::Metre;
  CHECK(read_one(to_step({box_shell(1.0)}, w), r).volume(kRel) == Approx(8.0).epsilon(1e-9));
  w.unit = StepLengthUnit::Inch;
  CHECK(read_one(to_step({box_shell(1.0)}, w)).volume(kRel) == Approx(8.0 * 25.4 * 25.4 * 25.4).epsilon(1e-9));
  r.unit = StepLengthUnit::Inch;
  CHECK(read_one(to_step({box_shell(1.0)}, w), r).volume(kRel) == Approx(8.0).epsilon(1e-9));
}

TEST_CASE("STEP read: cube from planes and lines", "[step_read]") {
  const Shell s = read_one(wrap(kCube));
  CHECK(s.faces().size() == 6);
  CHECK(s.edges().size() == 12);
  CHECK(s.is_closed());
  CHECK(s.volume(kRel) == Approx(1.0).epsilon(1e-12));
  CHECK(s.area() == Approx(6.0).epsilon(1e-10));
  // file in metres, read as millimetres
  const Shell m = read_one(wrap(kCube, "(LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT($,.METRE.))"));
  CHECK(m.volume(kRel) == Approx(1e9).epsilon(1e-12));
  // inch file (conversion based unit)
  const std::string inch =
      "(CONVERSION_BASED_UNIT('INCH',#902) LENGTH_UNIT() NAMED_UNIT(#903));\n"
      "#902=LENGTH_MEASURE_WITH_UNIT(LENGTH_MEASURE(25.4),#904);\n#903=DIMENSIONAL_EXPONENTS(1.,0.,0.,0.,0.,0.,0.);\n"
      "#904=(LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT(.MILLI.,.METRE.))";
  CHECK(read_one(wrap(kCube, inch)).volume(kRel) == Approx(25.4 * 25.4 * 25.4).epsilon(1e-12));
}

TEST_CASE("STEP read: cylinder from CYLINDRICAL_SURFACE and circles", "[step_read]") {
  const Shell s = read_one(wrap(kCylinder));
  CHECK(s.faces().size() == 3);
  CHECK(s.is_closed());
  CHECK(s.volume(kRel) == Approx(pi / 2.0).epsilon(1e-7));
  CHECK(s.area() == Approx(2 * pi * 0.5 * 2.0 + 2 * pi * 0.25).epsilon(1e-7));
}

TEST_CASE("STEP read: partial cylindrical face", "[step_read]") {
  const Shell s = read_one(wrap(kHalfCylinder));
  CHECK(s.faces().size() == 4);
  CHECK(s.is_closed());
  CHECK(s.volume(kRel) == Approx(pi / 2.0).epsilon(1e-7));
}

TEST_CASE("STEP read: spheres from SPHERICAL_SURFACE", "[step_read]") {
  const Shell s = read_one(wrap(kSphere));
  CHECK(s.faces().size() == 1);
  CHECK(s.is_closed());
  CHECK(s.volume(kRel) == Approx(4.0 / 3.0 * pi * 8.0).epsilon(1e-7));
  const Shell h = read_one(wrap(kHemisphere));
  CHECK(h.faces().size() == 2);
  CHECK(h.is_closed());
  CHECK(h.volume(kRel) == Approx(2.0 / 3.0 * pi).epsilon(1e-7));
  CHECK(h.area() == Approx(3.0 * pi).epsilon(1e-7));
}

TEST_CASE("STEP read: errors and skipped entities", "[step_read]") {
  CHECK_THROWS_AS(from_step("not a step file"), GeometryError);
  CHECK_THROWS_AS(from_step("ISO-10303-21;\nDATA;\n#1=CARTESIAN_POINT('',(0.,0.,0.);\nENDSEC;\n"), GeometryError);
  CHECK_THROWS_AS(read_step("does/not/exist.step"), std::runtime_error);
  CHECK(from_step(wrap("#1=CARTESIAN_POINT('',(0.,0.,0.));\n")).empty());
  // an unsupported surface type is named with its instance number
  std::string text = wrap(kCylinder);
  const auto pos = text.find("#30=CYLINDRICAL_SURFACE('',#7,0.5)");
  REQUIRE(pos != std::string::npos);
  text.replace(pos, 34, "#30=OFFSET_SURFACE('',#7,0.5,.F.)");
  try {
    (void)from_step(text);
    FAIL("expected an exception");
  } catch (const GeometryError& e) {
    const std::string what = e.what();
    CHECK(what.find("OFFSET_SURFACE") != std::string::npos);
    CHECK(what.find("#30") != std::string::npos);
  }
  std::vector<std::string> warnings;
  StepReadOptions opt;
  opt.skip_unsupported = true;
  opt.warnings = &warnings;
  const Shell s = read_one(text, opt);
  CHECK(s.faces().size() == 2);
  CHECK_FALSE(s.is_closed());
  REQUIRE(warnings.size() == 1);
  CHECK(warnings[0].find("OFFSET_SURFACE") != std::string::npos);
}
