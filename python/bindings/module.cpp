// Python bindings for the Nurbforge kernel (module nurbforge._core).
#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "nurbforge/nurbforge.hpp"

namespace py = pybind11;
using namespace nf;

namespace {

py::array_t<double> to_numpy(const std::vector<Vec3>& pts) {
  py::array_t<double> arr({static_cast<py::ssize_t>(pts.size()), py::ssize_t{3}});
  auto r = arr.mutable_unchecked<2>();
  for (py::ssize_t i = 0; i < static_cast<py::ssize_t>(pts.size()); ++i) {
    const auto& p = pts[static_cast<std::size_t>(i)];
    r(i, 0) = p.x; r(i, 1) = p.y; r(i, 2) = p.z;
  }
  return arr;
}

std::vector<Vec3> from_numpy(const py::array_t<double, py::array::c_style | py::array::forcecast>& a) {
  if (a.ndim() != 2 || a.shape(1) != 3) throw py::value_error("expected an array of shape (n, 3)");
  auto r = a.unchecked<2>();
  std::vector<Vec3> out(static_cast<std::size_t>(a.shape(0)));
  for (py::ssize_t i = 0; i < a.shape(0); ++i) out[static_cast<std::size_t>(i)] = {r(i, 0), r(i, 1), r(i, 2)};
  return out;
}

}  // namespace

PYBIND11_MODULE(_core, m) {
  m.doc() = "Nurbforge - headless NURBS CAD kernel (C++ core)";
  m.attr("__version__") = NURBFORGE_VERSION;

  py::register_exception<GeometryError>(m, "GeometryError", PyExc_ValueError);
  py::register_exception<NotImplementedError>(m, "NotImplementedError_", PyExc_NotImplementedError);

  // ------------------------------------------------------------------ Vec3
  py::class_<Vec3>(m, "Vec3")
      .def(py::init<>())
      .def(py::init<double, double, double>(), py::arg("x"), py::arg("y"), py::arg("z"))
      .def(py::init([](const std::array<double, 3>& a) { return Vec3(a); }))
      .def_readwrite("x", &Vec3::x)
      .def_readwrite("y", &Vec3::y)
      .def_readwrite("z", &Vec3::z)
      .def(py::self + py::self)
      .def(py::self - py::self)
      .def(py::self * double())
      .def(double() * py::self)
      .def(-py::self)
      .def("__iter__", [](const Vec3& v) { return py::iter(py::make_tuple(v.x, v.y, v.z)); })
      .def("__getitem__", [](const Vec3& v, int i) {
        if (i < 0 || i > 2) throw py::index_error();
        return i == 0 ? v.x : (i == 1 ? v.y : v.z);
      })
      .def("__len__", [](const Vec3&) { return 3; })
      .def("__repr__", [](const Vec3& v) {
        return "Vec3(" + std::to_string(v.x) + ", " + std::to_string(v.y) + ", " + std::to_string(v.z) + ")";
      })
      .def("to_tuple", [](const Vec3& v) { return py::make_tuple(v.x, v.y, v.z); });
  py::implicitly_convertible<py::tuple, Vec3>();
  py::implicitly_convertible<py::list, Vec3>();

  m.def("dot", &dot);
  m.def("cross", &cross);
  m.def("norm", &norm);
  m.def("distance", &distance);

  // ------------------------------------------------------------------ knots
  m.def("clamped_uniform_knots", &clamped_uniform_knots, py::arg("degree"), py::arg("n_ctrl"));
  m.def("unique_knots", &unique_knots, py::arg("knots"),
        "Distinct knot values and multiplicities (gmsh/OCC format).");
  m.def("validate_knots", &validate_knots, py::arg("knots"), py::arg("degree"), py::arg("n_ctrl"));

  // ------------------------------------------------------------------ tolerances
  py::class_<Tolerance>(m, "Tolerance", "Per-operation tolerances (length, parametric, angular).")
      .def(py::init([](double length, double parametric, double angular) {
             return Tolerance{length, parametric, angular};
           }),
           py::arg("length") = tol::kLength, py::arg("parametric") = tol::kParametric,
           py::arg("angular") = tol::kAngular)
      .def_readwrite("length", &Tolerance::length)
      .def_readwrite("parametric", &Tolerance::parametric)
      .def_readwrite("angular", &Tolerance::angular)
      .def("__repr__", [](const Tolerance& t) {
        return "Tolerance(length=" + std::to_string(t.length) + ", parametric=" +
               std::to_string(t.parametric) + ", angular=" + std::to_string(t.angular) + ")";
      });

  // ------------------------------------------------------------------ Transform
  py::class_<Transform>(m, "Transform")
      .def(py::init<>())
      .def(py::init<const std::array<double, 12>&>(), py::arg("matrix_3x4_row_major"))
      .def_static("translation", &Transform::translation, py::arg("d"))
      .def_static("rotation", &Transform::rotation, py::arg("origin"), py::arg("axis"), py::arg("angle"))
      .def_static("scaling", py::overload_cast<const Vec3&, double>(&Transform::scaling),
                  py::arg("center"), py::arg("factor"))
      .def_static("scaling", py::overload_cast<const Vec3&, const Vec3&>(&Transform::scaling),
                  py::arg("center"), py::arg("factors"))
      .def_static("mirror", &Transform::mirror, py::arg("point"), py::arg("normal"))
      .def("apply", &Transform::apply)
      .def("apply_vector", &Transform::apply_vector)
      .def("matrix", &Transform::matrix)
      .def(py::self * py::self);

  // ------------------------------------------------------------------ NurbsCurve
  py::class_<NurbsCurve>(m, "NurbsCurve")
      .def(py::init([](std::size_t degree, KnotVector knots,
                       const py::array_t<double, py::array::c_style | py::array::forcecast>& pts,
                       std::vector<double> weights) {
             return NurbsCurve(degree, std::move(knots), from_numpy(pts), std::move(weights));
           }),
           py::arg("degree"), py::arg("knots"), py::arg("control_points"),
           py::arg("weights") = std::vector<double>{})
      .def_static("from_points", [](std::size_t degree,
                                    const py::array_t<double, py::array::c_style | py::array::forcecast>& pts) {
             return NurbsCurve::from_points(degree, from_numpy(pts));
           }, py::arg("degree"), py::arg("control_points"))
      .def_property_readonly("degree", &NurbsCurve::degree)
      .def_property_readonly("knots", &NurbsCurve::knots)
      .def_property_readonly("weights", &NurbsCurve::weights)
      .def_property_readonly("control_points", [](const NurbsCurve& c) { return to_numpy(c.control_points()); })
      .def_property_readonly("is_rational", &NurbsCurve::is_rational)
      .def_property_readonly("domain", [](const NurbsCurve& c) {
        return py::make_tuple(c.domain_start(), c.domain_end());
      })
      .def("evaluate", &NurbsCurve::evaluate, py::arg("u"))
      .def("__call__", &NurbsCurve::evaluate, py::arg("u"))
      .def("derivatives", &NurbsCurve::derivatives, py::arg("u"), py::arg("order") = 1)
      .def("tangent", &NurbsCurve::tangent, py::arg("u"))
      .def("sample", [](const NurbsCurve& c, std::size_t n) { return to_numpy(c.sample(n)); },
           py::arg("n") = 100)
      .def("insert_knot", &NurbsCurve::insert_knot, py::arg("u"), py::arg("times") = 1)
      .def("refine_knots", &NurbsCurve::refine_knots, py::arg("knots"),
           "Insert all given knots at once (shape-preserving).")
      .def_property_readonly("is_clamped", &NurbsCurve::is_clamped)
      .def("clamped", &NurbsCurve::clamped, "Same curve with a clamped knot vector.")
      .def("elevate_degree", &NurbsCurve::elevate_degree, py::arg("times") = 1,
           "Degree elevation (shape-preserving); result is clamped.")
      .def("try_remove_knot", &NurbsCurve::try_remove_knot, py::arg("u"), py::arg("times") = 1,
           py::arg("tolerance") = tol::kLength,
           "Remove knot u up to `times` times within tolerance -> (curve, removed).")
      .def("remove_knot", &NurbsCurve::remove_knot, py::arg("u"), py::arg("times") = 1,
           py::arg("tolerance") = tol::kLength)
      .def("reduce_degree", &NurbsCurve::reduce_degree, py::arg("tolerance") = tol::kLength,
           "Degree reduction by one within tolerance (raises GeometryError otherwise).")
      .def("split", &NurbsCurve::split, py::arg("u"),
           "Split at u -> (left, right); both keep their parameter intervals.")
      .def("sub_curve", &NurbsCurve::sub_curve, py::arg("u0"), py::arg("u1"))
      .def("decompose_bezier", &NurbsCurve::decompose_bezier,
           "Bezier segments, each on its original parameter interval.")
      .def("reversed", &NurbsCurve::reversed)
      .def("transformed", &NurbsCurve::transformed, py::arg("transform"))
      .def("__repr__", [](const NurbsCurve& c) {
        return "<NurbsCurve degree=" + std::to_string(c.degree()) +
               " n_ctrl=" + std::to_string(c.num_control_points()) +
               (c.is_rational() ? " rational>" : ">");
      });

  // ------------------------------------------------------------------ NurbsSurface
  py::class_<NurbsSurface>(m, "NurbsSurface")
      .def(py::init([](std::size_t du, std::size_t dv, KnotVector ku, KnotVector kv, std::size_t nu,
                       std::size_t nv,
                       const py::array_t<double, py::array::c_style | py::array::forcecast>& pts,
                       std::vector<double> w) {
             return NurbsSurface(du, dv, std::move(ku), std::move(kv), nu, nv, from_numpy(pts), std::move(w));
           }),
           py::arg("degree_u"), py::arg("degree_v"), py::arg("knots_u"), py::arg("knots_v"),
           py::arg("num_u"), py::arg("num_v"), py::arg("control_points"),
           py::arg("weights") = std::vector<double>{},
           "Control points row-major: index = i * num_v + j (i along u).")
      .def_property_readonly("degree_u", &NurbsSurface::degree_u)
      .def_property_readonly("degree_v", &NurbsSurface::degree_v)
      .def_property_readonly("knots_u", &NurbsSurface::knots_u)
      .def_property_readonly("knots_v", &NurbsSurface::knots_v)
      .def_property_readonly("num_u", &NurbsSurface::num_u)
      .def_property_readonly("num_v", &NurbsSurface::num_v)
      .def_property_readonly("weights", &NurbsSurface::weights)
      .def_property_readonly("control_points", [](const NurbsSurface& s) { return to_numpy(s.control_points()); })
      .def_property_readonly("domain_u", &NurbsSurface::domain_u)
      .def_property_readonly("domain_v", &NurbsSurface::domain_v)
      .def("evaluate", &NurbsSurface::evaluate, py::arg("u"), py::arg("v"))
      .def("__call__", &NurbsSurface::evaluate, py::arg("u"), py::arg("v"))
      .def("derivatives", &NurbsSurface::derivatives, py::arg("u"), py::arg("v"), py::arg("order") = 1,
           "Partial derivatives D[k][l] = d^(k+l)S/du^k dv^l for k + l <= order.")
      .def("first_derivatives", &NurbsSurface::first_derivatives, py::arg("u"), py::arg("v"))
      .def("normal", &NurbsSurface::normal, py::arg("u"), py::arg("v"))
      .def("iso_curve_u", &NurbsSurface::iso_curve_u, py::arg("u"))
      .def("iso_curve_v", &NurbsSurface::iso_curve_v, py::arg("v"))
      .def("sample", [](const NurbsSurface& s, std::size_t nu, std::size_t nv) {
             auto a = to_numpy(s.sample(nu, nv));
             return a.reshape({static_cast<py::ssize_t>(nu), static_cast<py::ssize_t>(nv), py::ssize_t{3}});
           }, py::arg("nu") = 20, py::arg("nv") = 20)
      .def("refine_knots_u", &NurbsSurface::refine_knots_u, py::arg("knots"))
      .def("refine_knots_v", &NurbsSurface::refine_knots_v, py::arg("knots"))
      .def("elevate_degree_u", &NurbsSurface::elevate_degree_u, py::arg("times") = 1)
      .def("elevate_degree_v", &NurbsSurface::elevate_degree_v, py::arg("times") = 1)
      .def("clamped", &NurbsSurface::clamped)
      .def("split_u", &NurbsSurface::split_u, py::arg("u"))
      .def("split_v", &NurbsSurface::split_v, py::arg("v"))
      .def("sub_surface", &NurbsSurface::sub_surface, py::arg("u0"), py::arg("u1"), py::arg("v0"),
           py::arg("v1"))
      .def("try_remove_knot_u", &NurbsSurface::try_remove_knot_u, py::arg("u"), py::arg("times") = 1,
           py::arg("tolerance") = tol::kLength)
      .def("try_remove_knot_v", &NurbsSurface::try_remove_knot_v, py::arg("v"), py::arg("times") = 1,
           py::arg("tolerance") = tol::kLength)
      .def("remove_knot_u", &NurbsSurface::remove_knot_u, py::arg("u"), py::arg("times") = 1,
           py::arg("tolerance") = tol::kLength)
      .def("remove_knot_v", &NurbsSurface::remove_knot_v, py::arg("v"), py::arg("times") = 1,
           py::arg("tolerance") = tol::kLength)
      .def("reduce_degree_u", &NurbsSurface::reduce_degree_u, py::arg("tolerance") = tol::kLength)
      .def("reduce_degree_v", &NurbsSurface::reduce_degree_v, py::arg("tolerance") = tol::kLength)
      .def("decompose_bezier", &NurbsSurface::decompose_bezier,
           "Bezier patches in row-major order (u span outer), original parameter rectangles.")
      .def("transformed", &NurbsSurface::transformed, py::arg("transform"))
      .def("transposed", &NurbsSurface::transposed)
      .def("__repr__", [](const NurbsSurface& s) {
        return "<NurbsSurface degree=(" + std::to_string(s.degree_u()) + "," + std::to_string(s.degree_v()) +
               ") net=" + std::to_string(s.num_u()) + "x" + std::to_string(s.num_v()) + ">";
      });

  // ------------------------------------------------------------------ closest point
  py::class_<CurvePoint>(m, "CurvePoint")
      .def_readonly("u", &CurvePoint::u)
      .def_readonly("point", &CurvePoint::point)
      .def_readonly("distance", &CurvePoint::distance)
      .def("__repr__", [](const CurvePoint& r) {
        return "CurvePoint(u=" + std::to_string(r.u) + ", distance=" + std::to_string(r.distance) + ")";
      });
  py::class_<SurfacePoint>(m, "SurfacePoint")
      .def_readonly("u", &SurfacePoint::u)
      .def_readonly("v", &SurfacePoint::v)
      .def_readonly("point", &SurfacePoint::point)
      .def_readonly("distance", &SurfacePoint::distance)
      .def("__repr__", [](const SurfacePoint& r) {
        return "SurfacePoint(u=" + std::to_string(r.u) + ", v=" + std::to_string(r.v) +
               ", distance=" + std::to_string(r.distance) + ")";
      });
  m.def("closest_point", py::overload_cast<const NurbsCurve&, const Vec3&, const Tolerance&>(&closest_point),
        py::arg("curve"), py::arg("point"), py::arg("tolerance") = Tolerance{},
        "Global closest point on a curve (point projection / inversion).");
  m.def("closest_point", py::overload_cast<const NurbsSurface&, const Vec3&, const Tolerance&>(&closest_point),
        py::arg("surface"), py::arg("point"), py::arg("tolerance") = Tolerance{},
        "Global closest point on a surface.");
  m.def("closest_point_local",
        py::overload_cast<const NurbsCurve&, const Vec3&, double, const Tolerance&>(&closest_point_local),
        py::arg("curve"), py::arg("point"), py::arg("u0"), py::arg("tolerance") = Tolerance{});
  m.def("closest_point_local",
        py::overload_cast<const NurbsSurface&, const Vec3&, double, double, const Tolerance&>(
            &closest_point_local),
        py::arg("surface"), py::arg("point"), py::arg("u0"), py::arg("v0"),
        py::arg("tolerance") = Tolerance{});

  // ------------------------------------------------------------------ properties
  m.def("curvature", py::overload_cast<const NurbsCurve&, double>(&curvature), py::arg("curve"),
        py::arg("u"), "Curvature |C' x C''| / |C'|^3.");
  m.def("torsion", &torsion, py::arg("curve"), py::arg("u"));
  m.def("arc_length", py::overload_cast<const NurbsCurve&, double, double, double>(&arc_length),
        py::arg("curve"), py::arg("u0"), py::arg("u1"), py::arg("rel_tol") = kIntegrationRelTol);
  m.def("length", py::overload_cast<const NurbsCurve&, double>(&arc_length), py::arg("curve"),
        py::arg("rel_tol") = kIntegrationRelTol, "Total arc length of a curve.");
  m.def("parameter_at_length", &parameter_at_length, py::arg("curve"), py::arg("s"),
        py::arg("rel_tol") = kIntegrationRelTol);
  py::class_<SurfaceCurvature>(m, "SurfaceCurvature")
      .def_readonly("gaussian", &SurfaceCurvature::gaussian)
      .def_readonly("mean", &SurfaceCurvature::mean)
      .def_readonly("k1", &SurfaceCurvature::k1)
      .def_readonly("k2", &SurfaceCurvature::k2)
      .def_readonly("normal", &SurfaceCurvature::normal)
      .def_readonly("dir1", &SurfaceCurvature::dir1)
      .def_readonly("dir2", &SurfaceCurvature::dir2)
      .def("__repr__", [](const SurfaceCurvature& k) {
        return "SurfaceCurvature(K=" + std::to_string(k.gaussian) + ", H=" + std::to_string(k.mean) +
               ", k1=" + std::to_string(k.k1) + ", k2=" + std::to_string(k.k2) + ")";
      });
  m.def("surface_curvature", py::overload_cast<const NurbsSurface&, double, double>(&curvature),
        py::arg("surface"), py::arg("u"), py::arg("v"),
        "Gaussian, mean and principal curvatures and directions.");
  m.def("area", py::overload_cast<const NurbsSurface&, double>(&area), py::arg("surface"),
        py::arg("rel_tol") = kIntegrationRelTol);
  m.def("area", py::overload_cast<const NurbsSurface&, double, double, double, double, double>(&area),
        py::arg("surface"), py::arg("u0"), py::arg("u1"), py::arg("v0"), py::arg("v1"),
        py::arg("rel_tol") = kIntegrationRelTol);

  // ------------------------------------------------------------------ joining
  m.def("join", &join, py::arg("a"), py::arg("b"), py::arg("tolerance") = tol::kLength,
        "Join curve b to the end of curve a (gap <= tolerance).");
  m.def("join_u", &join_u, py::arg("a"), py::arg("b"), py::arg("tolerance") = tol::kLength);
  m.def("join_v", &join_v, py::arg("a"), py::arg("b"), py::arg("tolerance") = tol::kLength);

  // ------------------------------------------------------------------ primitives
  m.def("make_line", &make_line, py::arg("p0"), py::arg("p1"));
  m.def("make_polyline", &make_polyline, py::arg("points"));
  m.def("make_arc", &make_arc, py::arg("center"), py::arg("x_axis"), py::arg("y_axis"),
        py::arg("radius"), py::arg("start_angle"), py::arg("end_angle"));
  m.def("make_circle", &make_circle, py::arg("center"), py::arg("normal"), py::arg("radius"));
  m.def("make_bilinear", &make_bilinear, py::arg("p00"), py::arg("p10"), py::arg("p01"), py::arg("p11"));
  m.def("make_rectangle", &make_rectangle, py::arg("origin"), py::arg("edge_u"), py::arg("edge_v"));

  // ------------------------------------------------------------------ construction
  m.def("extrude", &extrude, py::arg("profile"), py::arg("direction"));
  m.def("revolve", &revolve, py::arg("profile"), py::arg("axis_point"), py::arg("axis_dir"),
        py::arg("angle"));
  m.def("ruled", &ruled, py::arg("c0"), py::arg("c1"));
  m.def("reparametrize", &reparametrize, py::arg("curve"), py::arg("a"), py::arg("b"));

  // ------------------------------------------------------------------ patterns
  m.def("linear_pattern", &linear_pattern, py::arg("step"), py::arg("count"));
  m.def("grid_pattern", &grid_pattern, py::arg("step_a"), py::arg("count_a"), py::arg("step_b"),
        py::arg("count_b"));
  m.def("polar_pattern", &polar_pattern, py::arg("axis_point"), py::arg("axis_dir"),
        py::arg("count"), py::arg("total_angle"));
}
