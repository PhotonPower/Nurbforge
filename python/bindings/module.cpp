// Python bindings for the Nurbforge kernel (module nurbforge._core).
#include <optional>
#include <tuple>

#include <pybind11/functional.h>
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

  // ------------------------------------------------------------------ intersections
  py::class_<CurveCurveIntersection>(m, "CurveCurveIntersection")
      .def_readonly("s", &CurveCurveIntersection::s)
      .def_readonly("t", &CurveCurveIntersection::t)
      .def_readonly("point", &CurveCurveIntersection::point)
      .def_readonly("distance", &CurveCurveIntersection::distance)
      .def("__repr__", [](const CurveCurveIntersection& x) {
        return "CurveCurveIntersection(s=" + std::to_string(x.s) + ", t=" + std::to_string(x.t) + ")";
      });
  py::class_<CurveSurfaceIntersection>(m, "CurveSurfaceIntersection")
      .def_readonly("t", &CurveSurfaceIntersection::t)
      .def_readonly("u", &CurveSurfaceIntersection::u)
      .def_readonly("v", &CurveSurfaceIntersection::v)
      .def_readonly("point", &CurveSurfaceIntersection::point)
      .def_readonly("distance", &CurveSurfaceIntersection::distance)
      .def("__repr__", [](const CurveSurfaceIntersection& x) {
        return "CurveSurfaceIntersection(t=" + std::to_string(x.t) + ", u=" + std::to_string(x.u) +
               ", v=" + std::to_string(x.v) + ")";
      });
  m.def("intersect",
        py::overload_cast<const NurbsCurve&, const NurbsCurve&, const Tolerance&>(&intersect),
        py::arg("a"), py::arg("b"), py::arg("tolerance") = Tolerance{},
        "All intersection points of two curves, sorted by the parameter on a.");
  m.def("intersect",
        py::overload_cast<const NurbsCurve&, const NurbsSurface&, const Tolerance&>(&intersect),
        py::arg("curve"), py::arg("surface"), py::arg("tolerance") = Tolerance{},
        "All intersection points of a curve and a surface, sorted by the curve parameter.");

  py::class_<SurfaceSurfacePoint>(m, "SurfaceSurfacePoint")
      .def_readonly("point", &SurfaceSurfacePoint::point)
      .def_readonly("u1", &SurfaceSurfacePoint::u1)
      .def_readonly("v1", &SurfaceSurfacePoint::v1)
      .def_readonly("u2", &SurfaceSurfacePoint::u2)
      .def_readonly("v2", &SurfaceSurfacePoint::v2);
  py::class_<IntersectionCurve>(m, "IntersectionCurve")
      .def_readonly("points", &IntersectionCurve::points)
      .def_readonly("closed", &IntersectionCurve::closed)
      .def_readonly("curve", &IntersectionCurve::curve)
      .def_readonly("pcurve1", &IntersectionCurve::pcurve1)
      .def_readonly("pcurve2", &IntersectionCurve::pcurve2)
      .def("__repr__", [](const IntersectionCurve& c) {
        return "<IntersectionCurve points=" + std::to_string(c.points.size()) +
               (c.closed ? " closed>" : " open>");
      });
  py::class_<SurfaceIntersection>(m, "SurfaceIntersection")
      .def_readonly("curves", &SurfaceIntersection::curves)
      .def_readonly("points", &SurfaceIntersection::points);
  py::class_<SurfaceIntersectionOptions>(m, "SurfaceIntersectionOptions")
      .def(py::init([](Tolerance tol, double chord, double fit) {
             return SurfaceIntersectionOptions{tol, chord, fit};
           }),
           py::arg("tolerance") = Tolerance{}, py::arg("chord_tolerance") = 0.0,
           py::arg("fit_tolerance") = 0.0)
      .def_readwrite("tolerance", &SurfaceIntersectionOptions::tol)
      .def_readwrite("chord_tolerance", &SurfaceIntersectionOptions::chord_tolerance)
      .def_readwrite("fit_tolerance", &SurfaceIntersectionOptions::fit_tolerance);
  m.def("intersect",
        py::overload_cast<const NurbsSurface&, const NurbsSurface&, const SurfaceIntersectionOptions&>(&intersect),
        py::arg("a"), py::arg("b"), py::arg("options") = SurfaceIntersectionOptions{},
        "Surface/surface intersection: curves (3D + p-curves on both surfaces) and touching points.");

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

  // ------------------------------------------------------------------ fitting
  py::enum_<Parametrization>(m, "Parametrization")
      .value("Uniform", Parametrization::Uniform)
      .value("ChordLength", Parametrization::ChordLength)
      .value("Centripetal", Parametrization::Centripetal);
  using Arr = py::array_t<double, py::array::c_style | py::array::forcecast>;
  const auto grid = [](const Arr& a) {
    if (a.ndim() != 3 || a.shape(2) != 3) throw py::value_error("expected a grid of shape (nu, nv, 3)");
    const auto nu = static_cast<std::size_t>(a.shape(0)), nv = static_cast<std::size_t>(a.shape(1));
    auto r = a.unchecked<3>();
    std::vector<Vec3> Q(nu * nv);
    for (py::ssize_t i = 0; i < a.shape(0); ++i)
      for (py::ssize_t j = 0; j < a.shape(1); ++j)
        Q[static_cast<std::size_t>(i) * nv + static_cast<std::size_t>(j)] = {r(i, j, 0), r(i, j, 1), r(i, j, 2)};
    return std::tuple{Q, nu, nv};
  };
  m.def("data_parameters", [](const Arr& pts, Parametrization method) {
        return data_parameters(from_numpy(pts), method);
      }, py::arg("points"), py::arg("method") = Parametrization::Centripetal);
  m.def("interpolate_curve",
        [](const Arr& pts, std::size_t degree, Parametrization method,
           const std::optional<std::pair<Vec3, Vec3>>& end_derivatives) {
          return interpolate_curve(from_numpy(pts), degree, method, end_derivatives);
        },
        py::arg("points"), py::arg("degree") = 3, py::arg("method") = Parametrization::Centripetal,
        py::arg("end_derivatives") = py::none(),
        "Global interpolation through points (n, 3); optional (C'(0), C'(1)).");
  m.def("approximate_curve",
        [](const Arr& pts, std::size_t degree, std::size_t n_ctrl, Parametrization method) {
          return approximate_curve(from_numpy(pts), degree, n_ctrl, method);
        },
        py::arg("points"), py::arg("degree"), py::arg("n_ctrl"),
        py::arg("method") = Parametrization::Centripetal);
  m.def("approximate_curve_tol",
        [](const Arr& pts, std::size_t degree, double tolerance, Parametrization method) {
          return approximate_curve_tol(from_numpy(pts), degree, tolerance, method);
        },
        py::arg("points"), py::arg("degree"), py::arg("tolerance"),
        py::arg("method") = Parametrization::Centripetal,
        "Least squares with the fewest control points meeting the tolerance.");
  m.def("interpolate_surface",
        [grid](const Arr& pts, std::size_t du, std::size_t dv, Parametrization method) {
          auto [Q, nu, nv] = grid(pts);
          return interpolate_surface(Q, nu, nv, du, dv, method);
        },
        py::arg("points"), py::arg("degree_u") = 3, py::arg("degree_v") = 3,
        py::arg("method") = Parametrization::Centripetal, "Interpolate a point grid (nu, nv, 3).");
  m.def("approximate_surface",
        [grid](const Arr& pts, std::size_t du, std::size_t dv, std::size_t cu, std::size_t cv,
               Parametrization method) {
          auto [Q, nu, nv] = grid(pts);
          return approximate_surface(Q, nu, nv, du, dv, cu, cv, method);
        },
        py::arg("points"), py::arg("degree_u"), py::arg("degree_v"), py::arg("ctrl_u"),
        py::arg("ctrl_v"), py::arg("method") = Parametrization::Centripetal);

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

  // ------------------------------------------------------------------ extrapolation
  py::enum_<Extension>(m, "Extension")
      .value("Natural", Extension::Natural)
      .value("Linear", Extension::Linear);
  py::enum_<CurveEnd>(m, "CurveEnd").value("Start", CurveEnd::Start).value("End", CurveEnd::End);
  py::enum_<SurfaceSide>(m, "SurfaceSide")
      .value("UMin", SurfaceSide::UMin)
      .value("UMax", SurfaceSide::UMax)
      .value("VMin", SurfaceSide::VMin)
      .value("VMax", SurfaceSide::VMax);
  m.def("extend", py::overload_cast<const NurbsCurve&, CurveEnd, double, Extension>(&extend),
        py::arg("curve"), py::arg("end"), py::arg("length"), py::arg("type") = Extension::Natural,
        "Extend a curve at one end by an arc length; original parameters are kept.");
  m.def("extend", py::overload_cast<const NurbsSurface&, SurfaceSide, double, Extension>(&extend),
        py::arg("surface"), py::arg("side"), py::arg("length"), py::arg("type") = Extension::Natural,
        "Extend a surface across one side; original parameters are kept.");

  // ------------------------------------------------------------------ loft
  m.def("loft_parameters", &loft_parameters, py::arg("sections"),
        py::arg("method") = Parametrization::Centripetal);
  m.def("loft", py::overload_cast<const std::vector<NurbsCurve>&, std::size_t, Parametrization>(&loft),
        py::arg("sections"), py::arg("degree_v") = 3, py::arg("method") = Parametrization::Centripetal,
        "Skinned surface through section curves (u along the sections, v across them).");
  m.def("loft",
        py::overload_cast<const std::vector<NurbsCurve>&, std::size_t, const std::vector<double>&>(&loft),
        py::arg("sections"), py::arg("degree_v"), py::arg("v_params"));
  m.def("make_compatible", [](std::vector<NurbsCurve> curves) {
        make_compatible(curves);
        return curves;
      }, py::arg("curves"), "Common degree, knot vector and domain [0, 1] for all curves.");

  // ------------------------------------------------------------------ sweep
  py::enum_<SweepFrame>(m, "SweepFrame")
      .value("RotationMinimizing", SweepFrame::RotationMinimizing)
      .value("Fixed", SweepFrame::Fixed);
  py::class_<SweepOptions>(m, "SweepOptions")
      .def(py::init([](SweepFrame frame, double twist, double tolerance, std::size_t max_sections) {
             return SweepOptions{frame, twist, tolerance, max_sections};
           }),
           py::arg("frame") = SweepFrame::RotationMinimizing, py::arg("twist") = 0.0,
           py::arg("tolerance") = 1e-6, py::arg("max_sections") = 1025)
      .def_readwrite("frame", &SweepOptions::frame)
      .def_readwrite("twist", &SweepOptions::twist)
      .def_readwrite("tolerance", &SweepOptions::tolerance)
      .def_readwrite("max_sections", &SweepOptions::max_sections);
  m.def("sweep_translational", &sweep_translational, py::arg("profile"), py::arg("trajectory"),
        "Exact translational sweep S(u, v) = C(u) + T(v) - T(v0).");
  m.def("sweep", &sweep, py::arg("profile"), py::arg("trajectory"),
        py::arg("options") = SweepOptions{},
        "Sweep a profile (given at the trajectory start) along a trajectory.");
  m.def("sweep2", &sweep2, py::arg("profile"), py::arg("rail1"), py::arg("rail2"),
        py::arg("options") = SweepOptions{}, "Two-rail sweep with uniform scaling.");

  // ------------------------------------------------------------------ offsets
  m.def("offset_curve", &offset_curve, py::arg("curve"), py::arg("distance"),
        py::arg("plane_normal") = Vec3{0, 0, 1}, py::arg("tolerance") = 1e-7,
        "Planar offset c(t) + d * normalize(plane_normal x c'(t)) as a cubic B-spline fitted to "
        "the tolerance (parameters kept). Raises GeometryError at cusps (|d| >= radius of "
        "curvature on the concave side), corners or non-planar curves.");
  m.def("offset_surface", &offset_surface, py::arg("surface"), py::arg("distance"),
        py::arg("tolerance") = 1e-7,
        "Offset surface S(u,v) + d * N(u,v) as a quintic tensor-product spline fitted to the "
        "tolerance (parameters kept). Raises GeometryError if d reaches a radius of curvature.");

  // ------------------------------------------------------------------ Coons / Gordon
  m.def("coons", &coons, py::arg("c0"), py::arg("c1"), py::arg("d0"), py::arg("d1"),
        py::arg("tolerance") = 1e-7,
        "Coons patch: S(u,0)=c0, S(u,1)=c1, S(0,v)=d0, S(1,v)=d1 (orientations fixed automatically).");
  m.def("gordon", &gordon, py::arg("u_curves"), py::arg("v_curves"), py::arg("tolerance") = 1e-7,
        "Gordon surface through a consistent curve network.");

  // ------------------------------------------------------------------ optical surfaces
  py::class_<Asphere>(m, "Asphere", "Rotationally symmetric sag: conic + even (r^4, r^6, ...) + odd (r, r^2, ...) terms.")
      .def(py::init([](double curvature, double conic, std::vector<double> even, std::vector<double> odd) {
             return Asphere{curvature, conic, std::move(even), std::move(odd)};
           }),
           py::arg("curvature") = 0.0, py::arg("conic") = 0.0, py::arg("even") = std::vector<double>{},
           py::arg("odd") = std::vector<double>{})
      .def_readwrite("curvature", &Asphere::curvature)
      .def_readwrite("conic", &Asphere::conic)
      .def_readwrite("even", &Asphere::even)
      .def_readwrite("odd", &Asphere::odd)
      .def("sag", &Asphere::sag, py::arg("r"))
      .def("max_radius", &Asphere::max_radius)
      .def_property_readonly("is_conic", &Asphere::is_conic);
  py::enum_<ZernikeIndexing>(m, "ZernikeIndexing")
      .value("Fringe", ZernikeIndexing::Fringe)
      .value("Noll", ZernikeIndexing::Noll);
  m.def("zernike", &zernike, py::arg("j"), py::arg("indexing"), py::arg("rho"), py::arg("theta"));
  m.def("zernike_nm", [](std::size_t j, ZernikeIndexing ix) {
        int n = 0, mm = 0;
        zernike_nm(j, ix, n, mm);
        return py::make_tuple(n, mm);
      }, py::arg("j"), py::arg("indexing"));
  py::class_<ZernikeSag>(m, "ZernikeSag")
      .def(py::init([](Asphere base, double norm_radius, ZernikeIndexing ix, std::vector<double> c) {
             return ZernikeSag{std::move(base), norm_radius, ix, std::move(c)};
           }),
           py::arg("base") = Asphere{}, py::arg("norm_radius") = 1.0,
           py::arg("indexing") = ZernikeIndexing::Fringe, py::arg("coeffs") = std::vector<double>{})
      .def_readwrite("base", &ZernikeSag::base)
      .def_readwrite("norm_radius", &ZernikeSag::norm_radius)
      .def_readwrite("indexing", &ZernikeSag::indexing)
      .def_readwrite("coeffs", &ZernikeSag::coeffs)
      .def("sag", &ZernikeSag::sag, py::arg("x"), py::arg("y"));
  py::class_<XYPolynomialSag>(m, "XYPolynomialSag")
      .def(py::init([](Asphere base, double norm_radius, const std::vector<std::tuple<int, int, double>>& terms) {
             XYPolynomialSag s{std::move(base), norm_radius, {}};
             for (const auto& [i, j, c] : terms) s.terms.push_back({i, j, c});
             return s;
           }),
           py::arg("base") = Asphere{}, py::arg("norm_radius") = 1.0,
           py::arg("terms") = std::vector<std::tuple<int, int, double>>{},
           "terms: list of (i, j, c) for c * (x/norm)^i * (y/norm)^j")
      .def("sag", &XYPolynomialSag::sag, py::arg("x"), py::arg("y"));
  m.def("asphere_profile", &asphere_profile, py::arg("asphere"), py::arg("r_max"),
        py::arg("tolerance") = 1e-7);
  m.def("asphere_surface", &asphere_surface, py::arg("asphere"), py::arg("r_max"),
        py::arg("tolerance") = 1e-7, "Surface of revolution; exact for pure conics.");
  m.def("sag_surface", &sag_surface, py::arg("f"), py::arg("x0"), py::arg("x1"), py::arg("y0"),
        py::arg("y1"), py::arg("degree") = 3, py::arg("tolerance") = 1e-7, py::arg("max_spans") = 256,
        "Graph surface z = f(x, y) over a rectangle, fitted to the tolerance.");
  m.def("zernike_surface", &zernike_surface, py::arg("zernike"), py::arg("half_width") = 0.0,
        py::arg("tolerance") = 1e-7);
  m.def("xy_polynomial_surface", &xy_polynomial_surface, py::arg("poly"), py::arg("half_width") = 0.0,
        py::arg("tolerance") = 1e-7);

  // ------------------------------------------------------------------ patterns
  m.def("linear_pattern", &linear_pattern, py::arg("step"), py::arg("count"));
  m.def("grid_pattern", &grid_pattern, py::arg("step_a"), py::arg("count_a"), py::arg("step_b"),
        py::arg("count_b"));
  m.def("polar_pattern", &polar_pattern, py::arg("axis_point"), py::arg("axis_dir"),
        py::arg("count"), py::arg("total_angle"));
}
