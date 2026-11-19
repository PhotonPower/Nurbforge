// Python bindings for the Nurbforge kernel (module nurbforge._core).
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <optional>
#include <tuple>

#include "nurbforge/nurbforge.hpp"

namespace py = pybind11;
using namespace nf;

namespace {

py::array_t<double> to_numpy(const std::vector<Vec3>& pts) {
  py::array_t<double> arr({static_cast<py::ssize_t>(pts.size()), py::ssize_t{3}});
  auto r = arr.mutable_unchecked<2>();
  for (py::ssize_t i = 0; i < static_cast<py::ssize_t>(pts.size()); ++i) {
    const auto& p = pts[static_cast<std::size_t>(i)];
    r(i, 0) = p.x;
    r(i, 1) = p.y;
    r(i, 2) = p.z;
  }
  return arr;
}

std::vector<Vec3> from_numpy(
    const py::array_t<double, py::array::c_style | py::array::forcecast>& a) {
  if (a.ndim() != 2 || a.shape(1) != 3) throw py::value_error("expected an array of shape (n, 3)");
  auto r = a.unchecked<2>();
  std::vector<Vec3> out(static_cast<std::size_t>(a.shape(0)));
  for (py::ssize_t i = 0; i < a.shape(0); ++i)
    out[static_cast<std::size_t>(i)] = {r(i, 0), r(i, 1), r(i, 2)};
  return out;
}

}  // namespace

namespace {

// A planar region from Python: a PlanarSection, or a sequence whose items are SectionLoops,
// closed curves, or sequences of curves forming one loop.
nf::Region to_region(const py::object& obj) {
  using namespace nf;
  Region r;
  if (py::isinstance<PlanarSection>(obj)) {
    for (const auto& loop : obj.cast<const PlanarSection&>().loops) {
      std::vector<NurbsCurve> l;
      for (const auto& sc : loop.curves) l.push_back(sc.curve);
      r.push_back(std::move(l));
    }
    return r;
  }
  for (const auto& item : obj) {
    if (py::isinstance<SectionLoop>(item)) {
      std::vector<NurbsCurve> l;
      for (const auto& sc : item.cast<const SectionLoop&>().curves) l.push_back(sc.curve);
      r.push_back(std::move(l));
    } else if (py::isinstance<NurbsCurve>(item)) {
      r.push_back({item.cast<NurbsCurve>()});
    } else {
      r.push_back(item.cast<std::vector<NurbsCurve>>());
    }
  }
  return r;
}

}  // namespace

PYBIND11_MODULE(_core, m) {
  m.doc() = "Nurbforge - headless NURBS CAD kernel (C++ core)";
  m.attr("__version__") = NURBFORGE_VERSION;

  py::register_exception<GeometryError>(m, "GeometryError", PyExc_ValueError);
  py::register_exception<NotImplementedError>(m, "NotImplementedError_", PyExc_NotImplementedError);

  // ------------------------------------------------------------------ Vec3
  py::class_<Vec3>(m, "Vec3")
      .def(py::pickle([](const Vec3& v) { return py::make_tuple(v.x, v.y, v.z); },
                      [](const py::tuple& t) {
                        return Vec3{t[0].cast<double>(), t[1].cast<double>(), t[2].cast<double>()};
                      }))
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
      .def("__getitem__",
           [](const Vec3& v, int i) {
             if (i < 0 || i > 2) throw py::index_error();
             return i == 0 ? v.x : (i == 1 ? v.y : v.z);
           })
      .def("__len__", [](const Vec3&) { return 3; })
      .def("__repr__",
           [](const Vec3& v) {
             return "Vec3(" + std::to_string(v.x) + ", " + std::to_string(v.y) + ", " +
                    std::to_string(v.z) + ")";
           })
      .def("to_tuple", [](const Vec3& v) { return py::make_tuple(v.x, v.y, v.z); });
  py::implicitly_convertible<py::tuple, Vec3>();
  py::implicitly_convertible<py::list, Vec3>();
  // any other sequence of three numbers, e.g. a NumPy row (np.array([x, y, z]))
  py::implicitly_convertible<py::sequence, Vec3>();

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
      .def(py::pickle(
          [](const Tolerance& t) { return py::make_tuple(t.length, t.parametric, t.angular); },
          [](const py::tuple& t) {
            return Tolerance{t[0].cast<double>(), t[1].cast<double>(), t[2].cast<double>()};
          }))
      .def(py::init([](double length, double parametric, double angular) {
             return Tolerance{length, parametric, angular};
           }),
           py::arg("length") = tol::kLength, py::arg("parametric") = tol::kParametric,
           py::arg("angular") = tol::kAngular)
      .def_readwrite("length", &Tolerance::length)
      .def_readwrite("parametric", &Tolerance::parametric)
      .def_readwrite("angular", &Tolerance::angular)
      .def("__repr__", [](const Tolerance& t) {
        return "Tolerance(length=" + std::to_string(t.length) +
               ", parametric=" + std::to_string(t.parametric) +
               ", angular=" + std::to_string(t.angular) + ")";
      });

  // ------------------------------------------------------------------ Transform
  py::class_<Transform>(m, "Transform")
      .def(py::pickle(
          [](const Transform& t) { return py::make_tuple(t.matrix()); },
          [](const py::tuple& t) { return Transform(t[0].cast<std::array<double, 12>>()); }))
      .def(py::init<>())
      .def(py::init<const std::array<double, 12>&>(), py::arg("matrix_3x4_row_major"))
      .def_static("translation", &Transform::translation, py::arg("d"))
      .def_static("rotation", &Transform::rotation, py::arg("origin"), py::arg("axis"),
                  py::arg("angle"))
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
      .def(py::pickle(
          [](const NurbsCurve& c) {
            return py::make_tuple(c.degree(), c.knots(), c.control_points(), c.weights());
          },
          [](const py::tuple& t) {
            return NurbsCurve(t[0].cast<std::size_t>(), t[1].cast<KnotVector>(),
                              t[2].cast<std::vector<Vec3>>(), t[3].cast<std::vector<double>>());
          }))
      .def(py::init([](std::size_t degree, KnotVector knots,
                       const py::array_t<double, py::array::c_style | py::array::forcecast>& pts,
                       std::vector<double> weights) {
             return NurbsCurve(degree, std::move(knots), from_numpy(pts), std::move(weights));
           }),
           py::arg("degree"), py::arg("knots"), py::arg("control_points"),
           py::arg("weights") = std::vector<double>{})
      .def_static(
          "from_points",
          [](std::size_t degree,
             const py::array_t<double, py::array::c_style | py::array::forcecast>& pts) {
            return NurbsCurve::from_points(degree, from_numpy(pts));
          },
          py::arg("degree"), py::arg("control_points"))
      .def_property_readonly("degree", &NurbsCurve::degree)
      .def_property_readonly("knots", &NurbsCurve::knots)
      .def_property_readonly("weights", &NurbsCurve::weights)
      .def_property_readonly("control_points",
                             [](const NurbsCurve& c) { return to_numpy(c.control_points()); })
      .def_property_readonly("is_rational", &NurbsCurve::is_rational)
      .def_property_readonly(
          "domain",
          [](const NurbsCurve& c) { return py::make_tuple(c.domain_start(), c.domain_end()); })
      .def("evaluate", &NurbsCurve::evaluate, py::arg("u"))
      .def("__call__", &NurbsCurve::evaluate, py::arg("u"))
      .def("derivatives", &NurbsCurve::derivatives, py::arg("u"), py::arg("order") = 1)
      .def("tangent", &NurbsCurve::tangent, py::arg("u"))
      .def(
          "sample", [](const NurbsCurve& c, std::size_t n) { return to_numpy(c.sample(n)); },
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
      .def("__rmul__", [](const NurbsCurve& c, const Transform& t) { return c.transformed(t); })
      .def("__repr__", [](const NurbsCurve& c) {
        return "<NurbsCurve degree=" + std::to_string(c.degree()) +
               " n_ctrl=" + std::to_string(c.num_control_points()) +
               (c.is_rational() ? " rational>" : ">");
      });

  // ------------------------------------------------------------------ NurbsSurface
  py::class_<NurbsSurface>(m, "NurbsSurface")
      .def(py::pickle(
          [](const NurbsSurface& s) {
            return py::make_tuple(s.degree_u(), s.degree_v(), s.knots_u(), s.knots_v(), s.num_u(),
                                  s.num_v(), s.control_points(), s.weights());
          },
          [](const py::tuple& t) {
            return NurbsSurface(t[0].cast<std::size_t>(), t[1].cast<std::size_t>(),
                                t[2].cast<KnotVector>(), t[3].cast<KnotVector>(),
                                t[4].cast<std::size_t>(), t[5].cast<std::size_t>(),
                                t[6].cast<std::vector<Vec3>>(), t[7].cast<std::vector<double>>());
          }))
      .def(py::init([](std::size_t du, std::size_t dv, KnotVector ku, KnotVector kv, std::size_t nu,
                       std::size_t nv,
                       const py::array_t<double, py::array::c_style | py::array::forcecast>& pts,
                       std::vector<double> w) {
             return NurbsSurface(du, dv, std::move(ku), std::move(kv), nu, nv, from_numpy(pts),
                                 std::move(w));
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
      .def_property_readonly("control_points",
                             [](const NurbsSurface& s) { return to_numpy(s.control_points()); })
      .def_property_readonly("domain_u", &NurbsSurface::domain_u)
      .def_property_readonly("domain_v", &NurbsSurface::domain_v)
      .def("evaluate", &NurbsSurface::evaluate, py::arg("u"), py::arg("v"))
      .def("__call__", &NurbsSurface::evaluate, py::arg("u"), py::arg("v"))
      .def("derivatives", &NurbsSurface::derivatives, py::arg("u"), py::arg("v"),
           py::arg("order") = 1,
           "Partial derivatives D[k][l] = d^(k+l)S/du^k dv^l for k + l <= order.")
      .def("first_derivatives", &NurbsSurface::first_derivatives, py::arg("u"), py::arg("v"))
      .def("normal", &NurbsSurface::normal, py::arg("u"), py::arg("v"))
      .def("iso_curve_u", &NurbsSurface::iso_curve_u, py::arg("u"))
      .def("iso_curve_v", &NurbsSurface::iso_curve_v, py::arg("v"))
      .def(
          "sample",
          [](const NurbsSurface& s, std::size_t nu, std::size_t nv) {
            auto a = to_numpy(s.sample(nu, nv));
            return a.reshape(
                {static_cast<py::ssize_t>(nu), static_cast<py::ssize_t>(nv), py::ssize_t{3}});
          },
          py::arg("nu") = 20, py::arg("nv") = 20)
      .def("refine_knots_u", &NurbsSurface::refine_knots_u, py::arg("knots"))
      .def("refine_knots_v", &NurbsSurface::refine_knots_v, py::arg("knots"))
      .def("elevate_degree_u", &NurbsSurface::elevate_degree_u, py::arg("times") = 1)
      .def("elevate_degree_v", &NurbsSurface::elevate_degree_v, py::arg("times") = 1)
      .def("clamped", &NurbsSurface::clamped)
      .def("split_u", &NurbsSurface::split_u, py::arg("u"))
      .def("split_v", &NurbsSurface::split_v, py::arg("v"))
      .def("sub_surface", &NurbsSurface::sub_surface, py::arg("u0"), py::arg("u1"), py::arg("v0"),
           py::arg("v1"))
      .def("try_remove_knot_u", &NurbsSurface::try_remove_knot_u, py::arg("u"),
           py::arg("times") = 1, py::arg("tolerance") = tol::kLength)
      .def("try_remove_knot_v", &NurbsSurface::try_remove_knot_v, py::arg("v"),
           py::arg("times") = 1, py::arg("tolerance") = tol::kLength)
      .def("remove_knot_u", &NurbsSurface::remove_knot_u, py::arg("u"), py::arg("times") = 1,
           py::arg("tolerance") = tol::kLength)
      .def("remove_knot_v", &NurbsSurface::remove_knot_v, py::arg("v"), py::arg("times") = 1,
           py::arg("tolerance") = tol::kLength)
      .def("reduce_degree_u", &NurbsSurface::reduce_degree_u, py::arg("tolerance") = tol::kLength)
      .def("reduce_degree_v", &NurbsSurface::reduce_degree_v, py::arg("tolerance") = tol::kLength)
      .def("decompose_bezier", &NurbsSurface::decompose_bezier,
           "Bezier patches in row-major order (u span outer), original parameter rectangles.")
      .def("transformed", &NurbsSurface::transformed, py::arg("transform"))
      .def("__rmul__", [](const NurbsSurface& s, const Transform& t) { return s.transformed(t); })
      .def("transposed", &NurbsSurface::transposed)
      .def("__repr__", [](const NurbsSurface& s) {
        return "<NurbsSurface degree=(" + std::to_string(s.degree_u()) + "," +
               std::to_string(s.degree_v()) + ") net=" + std::to_string(s.num_u()) + "x" +
               std::to_string(s.num_v()) + ">";
      });

  // ------------------------------------------------------------------ closest point
  py::class_<CurvePoint>(m, "CurvePoint")
      .def_readonly("u", &CurvePoint::u)
      .def_readonly("point", &CurvePoint::point)
      .def_readonly("distance", &CurvePoint::distance)
      .def("__repr__", [](const CurvePoint& r) {
        return "CurvePoint(u=" + std::to_string(r.u) + ", distance=" + std::to_string(r.distance) +
               ")";
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
  m.def("closest_point",
        py::overload_cast<const NurbsCurve&, const Vec3&, const Tolerance&>(&closest_point),
        py::call_guard<py::gil_scoped_release>(), py::arg("curve"), py::arg("point"),
        py::arg("tolerance") = Tolerance{},
        "Global closest point on a curve (point projection / inversion).");
  m.def("closest_point",
        py::overload_cast<const NurbsSurface&, const Vec3&, const Tolerance&>(&closest_point),
        py::call_guard<py::gil_scoped_release>(), py::arg("surface"), py::arg("point"),
        py::arg("tolerance") = Tolerance{}, "Global closest point on a surface.");
  m.def("closest_point_local",
        py::overload_cast<const NurbsCurve&, const Vec3&, double, const Tolerance&>(
            &closest_point_local),
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
        return "CurveCurveIntersection(s=" + std::to_string(x.s) + ", t=" + std::to_string(x.t) +
               ")";
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
        py::call_guard<py::gil_scoped_release>(), py::arg("a"), py::arg("b"),
        py::arg("tolerance") = Tolerance{},
        "All intersection points of two curves, sorted by the parameter on a.");
  m.def("intersect",
        py::overload_cast<const NurbsCurve&, const NurbsSurface&, const Tolerance&>(&intersect),
        py::call_guard<py::gil_scoped_release>(), py::arg("curve"), py::arg("surface"),
        py::arg("tolerance") = Tolerance{},
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
      .def(py::pickle(
          [](const SurfaceIntersectionOptions& o) {
            return py::make_tuple(o.tol, o.chord_tolerance, o.fit_tolerance);
          },
          [](const py::tuple& t) {
            SurfaceIntersectionOptions o;
            o.tol = t[0].cast<Tolerance>();
            o.chord_tolerance = t[1].cast<double>();
            o.fit_tolerance = t[2].cast<double>();
            return o;
          }))
      .def(py::init([](Tolerance tol, double chord, double fit) {
             return SurfaceIntersectionOptions{tol, chord, fit};
           }),
           py::arg("tolerance") = Tolerance{}, py::arg("chord_tolerance") = 0.0,
           py::arg("fit_tolerance") = 0.0)
      .def_readwrite("tolerance", &SurfaceIntersectionOptions::tol)
      .def_readwrite("chord_tolerance", &SurfaceIntersectionOptions::chord_tolerance)
      .def_readwrite("fit_tolerance", &SurfaceIntersectionOptions::fit_tolerance);
  m.def(
      "intersect",
      py::overload_cast<const NurbsSurface&, const NurbsSurface&,
                        const SurfaceIntersectionOptions&>(&intersect),
      py::call_guard<py::gil_scoped_release>(), py::arg("a"), py::arg("b"),
      py::arg("options") = SurfaceIntersectionOptions{},
      "Surface/surface intersection: curves (3D + p-curves on both surfaces) and touching points.");

  py::class_<Plane>(m, "Plane", "Infinite plane through `point` with normal `normal`.")
      .def(py::pickle(
          [](const Plane& p) { return py::make_tuple(p.point, p.normal); },
          [](const py::tuple& t) { return Plane{t[0].cast<Vec3>(), t[1].cast<Vec3>()}; }))
      .def(py::init([](const Vec3& point, const Vec3& normal) { return Plane{point, normal}; }),
           py::arg("point"), py::arg("normal") = Vec3{0, 0, 1})
      .def_readwrite("point", &Plane::point)
      .def_readwrite("normal", &Plane::normal)
      .def("distance", &Plane::distance, py::arg("p"));
  m.def("planar_section", &planar_section, py::call_guard<py::gil_scoped_release>(),
        py::arg("surface"), py::arg("plane"), py::arg("options") = SurfaceIntersectionOptions{},
        "Section of a surface with a plane (pcurve2 in the plane's in-plane coordinates).");
  m.def("planar_sections", &planar_sections, py::call_guard<py::gil_scoped_release>(),
        py::arg("surface"), py::arg("plane"), py::arg("spacing"), py::arg("count"),
        py::arg("options") = SurfaceIntersectionOptions{});
  // shell sections (slicing of B-rep bodies)
  py::class_<SectionOptions>(m, "SectionOptions")
      .def(py::pickle(
          [](const SectionOptions& o) {
            return py::make_tuple(o.join_tolerance, o.fit_tolerance, o.s_axis);
          },
          [](const py::tuple& t) {
            SectionOptions o;
            o.join_tolerance = t[0].cast<double>();
            o.fit_tolerance = t[1].cast<double>();
            o.s_axis = t[2].cast<Vec3>();
            return o;
          }))
      .def(py::init<>())
      .def_readwrite("join_tolerance", &SectionOptions::join_tolerance)
      .def_readwrite("fit_tolerance", &SectionOptions::fit_tolerance)
      .def_readwrite("s_axis", &SectionOptions::s_axis);
  py::class_<SectionCurve>(m, "SectionCurve", "Piece of a section loop on one face of a shell.")
      .def_readonly("curve", &SectionCurve::curve, "In plane coordinates (s, t, 0).")
      .def_readonly("curve3d", &SectionCurve::curve3d)
      .def_readonly("pcurve", &SectionCurve::pcurve, "(u, v, 0) on the source face's surface.")
      .def_readonly("face", &SectionCurve::face)
      .def_readonly("exact", &SectionCurve::exact)
      .def("normal", &SectionCurve::normal, py::arg("t"),
           "Outward unit normal of the shell at parameter t.");
  py::class_<SectionLoop>(m, "SectionLoop")
      .def_readonly("curves", &SectionLoop::curves)
      .def_readonly("area", &SectionLoop::area, "Signed area in (s, t): > 0 outer loop, < 0 hole.");
  py::class_<PlanarSection>(m, "PlanarSection")
      .def_readonly("origin", &PlanarSection::origin)
      .def_readonly("s_axis", &PlanarSection::s_axis)
      .def_readonly("t_axis", &PlanarSection::t_axis)
      .def_readonly("normal", &PlanarSection::normal)
      .def_readonly("loops", &PlanarSection::loops);
  m.def("section", &section, py::call_guard<py::gil_scoped_release>(), py::arg("shell"),
        py::arg("plane"), py::arg("options") = SectionOptions{},
        "Section of a closed shell with a plane: closed loops oriented with the solid on the left "
        "(outer counter-clockwise, holes clockwise in the plane coordinates).");
  // hatching of section loops (contour and hatch lines of a layer)
  py::class_<HatchSegment>(m, "HatchSegment",
                           "Hatch segment inside the region bounded by the loops.")
      .def_readonly("start", &HatchSegment::start)
      .def_readonly("end", &HatchSegment::end)
      .def_readonly("line", &HatchSegment::line, "Index k of the line n . p = phase + k * spacing.")
      .def_readonly("region", &HatchSegment::region,
                    "Outer loop of the region containing the segment.")
      .def_readonly("start_loop", &HatchSegment::start_loop)
      .def_readonly("end_loop", &HatchSegment::end_loop)
      .def_property_readonly("length", &HatchSegment::length)
      .def("__repr__", [](const HatchSegment& s) {
        return "HatchSegment(line=" + std::to_string(s.line) +
               ", length=" + std::to_string(s.length()) + ", region=" + std::to_string(s.region) +
               ")";
      });
  py::class_<HatchContour>(m, "HatchContour",
                           "Closed contour: boundary of the region shrunk by `distance`.")
      .def_readonly("curves", &HatchContour::curves)
      .def_readonly("loop", &HatchContour::loop, "First of `loops`.")
      .def_readonly("loops", &HatchContour::loops,
                    "Loops it was offset from (merged contours: several).")
      .def_readonly("index", &HatchContour::index, "Index of the distance in `contours`.")
      .def_readonly("distance", &HatchContour::distance);
  m.def(
      "hatch",
      [](const py::object& loops, const py::sequence& direction, double spacing, double phase,
         const std::vector<double>& contours, double min_length, double tolerance,
         double contour_tolerance) {
        if (direction.size() != 2 && direction.size() != 3)
          throw py::value_error("direction must have 2 or 3 components");
        const Vec3 d{direction[0].cast<double>(), direction[1].cast<double>(), 0.0};
        HatchOptions o;
        o.phase = phase;
        o.contours = contours;
        o.min_length = min_length;
        o.tolerance = tolerance;
        o.contour_tolerance = contour_tolerance;
        Hatching h;
        if (py::isinstance<PlanarSection>(loops)) {
          h = hatch(loops.cast<const PlanarSection&>().loops, d, spacing, o);
        } else {
          std::vector<std::vector<NurbsCurve>> ls;
          for (const auto& item : loops) {
            if (py::isinstance<SectionLoop>(item)) {
              std::vector<NurbsCurve> l;
              for (const auto& sc : item.cast<const SectionLoop&>().curves) l.push_back(sc.curve);
              ls.push_back(std::move(l));
            } else if (py::isinstance<NurbsCurve>(item)) {
              ls.push_back({item.cast<NurbsCurve>()});
            } else {
              ls.push_back(item.cast<std::vector<NurbsCurve>>());
            }
          }
          py::gil_scoped_release release;
          h = hatch(ls, d, spacing, o);
        }
        return py::make_tuple(h.segments, h.contours);
      },
      py::arg("loops"), py::arg("direction"), py::arg("spacing"), py::kw_only(),
      py::arg("phase") = 0.0, py::arg("contours") = std::vector<double>{},
      py::arg("min_length") = 0.0, py::arg("tolerance") = 0.0, py::arg("contour_tolerance") = 0.0,
      "Hatch and contour lines of the region bounded by closed planar loops (region on the "
      "left).\n\n"
      "`loops`: a PlanarSection, its loops, or closed chains of curves in the xy plane (a single "
      "closed\n"
      "curve counts as a loop). Lines n . p = phase + k * spacing along `direction` (2 or 3\n"
      "components, z ignored), n = (-d.y, d.x). Returns (segments, contours): HatchSegment list "
      "sorted\n"
      "by line and position, and one HatchContour per loop and inward offset distance in "
      "`contours`.");
  m.def(
      "critical_heights", &critical_heights, py::arg("shell"), py::arg("direction"),
      py::arg("tolerance") = Tolerance{},
      "Sorted heights along `direction` at which the topology of the planar sections can change.");
  m.def("is_extrusion_band", &is_extrusion_band, py::call_guard<py::gil_scoped_release>(),
        py::arg("shell"), py::arg("direction"), py::arg("h0"), py::arg("h1"),
        py::arg("tolerance") = Tolerance{},
        "True if all faces meeting the band h0 < p . direction < h1 are parallel to `direction` "
        "there.");
  m.def("intersect",
        py::overload_cast<const NurbsCurve&, const Plane&, const Tolerance&>(&intersect),
        py::call_guard<py::gil_scoped_release>(), py::arg("curve"), py::arg("plane"),
        py::arg("tolerance") = Tolerance{}, "Intersection points of a curve with a plane.");

  // ------------------------------------------------------------------ properties
  m.def("curvature", py::overload_cast<const NurbsCurve&, double>(&curvature), py::arg("curve"),
        py::arg("u"), "Curvature ``|C' x C''| / |C'|^3``.");
  m.def("torsion", &torsion, py::arg("curve"), py::arg("u"));
  m.def("arc_length", py::overload_cast<const NurbsCurve&, double, double, double>(&arc_length),
        py::call_guard<py::gil_scoped_release>(), py::arg("curve"), py::arg("u0"), py::arg("u1"),
        py::arg("rel_tol") = kIntegrationRelTol);
  m.def("length", py::overload_cast<const NurbsCurve&, double>(&arc_length),
        py::call_guard<py::gil_scoped_release>(), py::arg("curve"),
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
        return "SurfaceCurvature(K=" + std::to_string(k.gaussian) +
               ", H=" + std::to_string(k.mean) + ", k1=" + std::to_string(k.k1) +
               ", k2=" + std::to_string(k.k2) + ")";
      });
  m.def("surface_curvature", py::overload_cast<const NurbsSurface&, double, double>(&curvature),
        py::arg("surface"), py::arg("u"), py::arg("v"),
        "Gaussian, mean and principal curvatures and directions.");
  m.def("area", py::overload_cast<const NurbsSurface&, double>(&area),
        py::call_guard<py::gil_scoped_release>(), py::arg("surface"),
        py::arg("rel_tol") = kIntegrationRelTol);
  m.def("area",
        py::overload_cast<const NurbsSurface&, double, double, double, double, double>(&area),
        py::call_guard<py::gil_scoped_release>(), py::arg("surface"), py::arg("u0"), py::arg("u1"),
        py::arg("v0"), py::arg("v1"), py::arg("rel_tol") = kIntegrationRelTol);

  // ------------------------------------------------------------------ fitting
  py::enum_<Parametrization>(m, "Parametrization")
      .value("Uniform", Parametrization::Uniform)
      .value("ChordLength", Parametrization::ChordLength)
      .value("Centripetal", Parametrization::Centripetal);
  using Arr = py::array_t<double, py::array::c_style | py::array::forcecast>;
  const auto grid = [](const Arr& a) {
    if (a.ndim() != 3 || a.shape(2) != 3)
      throw py::value_error("expected a grid of shape (nu, nv, 3)");
    const auto nu = static_cast<std::size_t>(a.shape(0)), nv = static_cast<std::size_t>(a.shape(1));
    auto r = a.unchecked<3>();
    std::vector<Vec3> Q(nu * nv);
    for (py::ssize_t i = 0; i < a.shape(0); ++i)
      for (py::ssize_t j = 0; j < a.shape(1); ++j)
        Q[static_cast<std::size_t>(i) * nv + static_cast<std::size_t>(j)] = {r(i, j, 0), r(i, j, 1),
                                                                             r(i, j, 2)};
    return std::tuple{Q, nu, nv};
  };
  m.def(
      "data_parameters",
      [](const Arr& pts, Parametrization method) {
        return data_parameters(from_numpy(pts), method);
      },
      py::arg("points"), py::arg("method") = Parametrization::Centripetal);
  m.def(
      "interpolate_curve",
      [](const Arr& pts, std::size_t degree, Parametrization method,
         const std::optional<std::pair<Vec3, Vec3>>& end_derivatives) {
        return interpolate_curve(from_numpy(pts), degree, method, end_derivatives);
      },
      py::arg("points"), py::arg("degree") = 3, py::arg("method") = Parametrization::Centripetal,
      py::arg("end_derivatives") = py::none(),
      "Global interpolation through points (n, 3); optional (C'(0), C'(1)).");
  py::class_<LineFit>(m, "LineFit", "Best-fit line (orthogonal regression).")
      .def_readonly("point", &LineFit::point, "centroid of the points")
      .def_readonly("direction", &LineFit::direction, "unit direction")
      .def_readonly("max_deviation", &LineFit::max_deviation,
                    "largest perpendicular distance of a point from the line")
      .def_readonly("rms", &LineFit::rms, "root mean square of the perpendicular distances")
      .def("__repr__", [](const LineFit& f) {
        return "<LineFit max_deviation=" + std::to_string(f.max_deviation) + ">";
      });
  py::class_<PlaneFit>(m, "PlaneFit", "Best-fit plane (orthogonal regression).")
      .def_readonly("plane", &PlaneFit::plane, "plane through the centroid, unit normal")
      .def_readonly("max_deviation", &PlaneFit::max_deviation,
                    "largest distance of a point from the plane")
      .def_readonly("rms", &PlaneFit::rms, "root mean square of the distances")
      .def("__repr__", [](const PlaneFit& f) {
        return "<PlaneFit max_deviation=" + std::to_string(f.max_deviation) + ">";
      });
  m.def(
      "fit_line", [](const Arr& pts) { return fit_line(from_numpy(pts)); }, py::arg("points"),
      "Best-fit line through points (n, 3) by orthogonal regression (minimal sum of squared "
      "perpendicular distances): centroid, unit direction, max and rms deviation.");
  m.def(
      "fit_plane", [](const Arr& pts) { return fit_plane(from_numpy(pts)); }, py::arg("points"),
      "Best-fit plane through points (n, 3) by orthogonal regression: plane through the "
      "centroid with unit normal, max and rms deviation.");
  py::class_<FittedSurface>(m, "FittedSurface", "Result of fit_surface.")
      .def_readonly("surface", &FittedSurface::surface, "the fitted surface")
      .def_readonly("max_error", &FittedSurface::max_error,
                    "largest distance of a point from the surface")
      .def_readonly("rms_error", &FittedSurface::rms_error, "rms of the point distances")
      .def_property_readonly(
          "parameters",
          [](const FittedSurface& f) {
            py::array_t<double> a({f.parameters.size(), std::size_t{2}});
            auto r = a.mutable_unchecked<2>();
            for (std::size_t k = 0; k < f.parameters.size(); ++k) {
              r(k, 0) = f.parameters[k][0];
              r(k, 1) = f.parameters[k][1];
            }
            return a;
          },
          "(n, 2) array: (u, v) of each point's foot point")
      .def("__repr__", [](const FittedSurface& f) {
        return "<FittedSurface max_error=" + std::to_string(f.max_error) +
               " rms_error=" + std::to_string(f.rms_error) + ">";
      });
  m.def(
      "fit_surface",
      [](const Arr& pts, std::size_t num_u, std::size_t num_v, std::size_t degree_u,
         std::size_t degree_v, double smoothing, std::size_t iterations,
         const std::optional<NurbsSurface>& base, const std::vector<NurbsCurve>& boundary) {
        SurfaceFitOptions o;
        o.num_u = num_u;
        o.num_v = num_v;
        o.degree_u = degree_u;
        o.degree_v = degree_v;
        o.smoothing = smoothing;
        o.iterations = iterations;
        o.base = base;
        o.boundary = boundary;
        const auto P = from_numpy(pts);
        py::gil_scoped_release release;
        return fit_surface(P, o);
      },
      py::arg("points"), py::arg("num_u") = 0, py::arg("num_v") = 0, py::arg("degree_u") = 3,
      py::arg("degree_v") = 3, py::arg("smoothing") = 1e-3, py::arg("iterations") = 2,
      py::arg("base") = std::nullopt, py::arg("boundary") = std::vector<NurbsCurve>{},
      "Fit a B-spline surface (domain [0, 1]^2) to unordered points (n, 3): least squares plus "
      "smoothing x thin-plate energy (relative weight; 1e-3 barely changes smooth data, 0.1-1 "
      "smooths noise), then `iterations` rounds of parameter correction. Parameters from the "
      "best-fit plane (height fields), from `base` (closest points), or from the Coons patch "
      "of four `boundary` curves [c0, c1, d0, d1], which are then interpolated exactly. "
      "num_u/num_v = 0: about sqrt(n / 4). Returns FittedSurface (surface, max_error, "
      "rms_error, parameters).");
  m.def(
      "approximate_curve",
      [](const Arr& pts, std::size_t degree, std::size_t n_ctrl, Parametrization method) {
        return approximate_curve(from_numpy(pts), degree, n_ctrl, method);
      },
      py::arg("points"), py::arg("degree"), py::arg("n_ctrl"),
      py::arg("method") = Parametrization::Centripetal);
  m.def(
      "approximate_curve_tol",
      [](const Arr& pts, std::size_t degree, double tolerance, Parametrization method) {
        return approximate_curve_tol(from_numpy(pts), degree, tolerance, method);
      },
      py::arg("points"), py::arg("degree"), py::arg("tolerance"),
      py::arg("method") = Parametrization::Centripetal,
      "Least squares with the fewest control points meeting the tolerance.");
  m.def(
      "interpolate_surface",
      [grid](const Arr& pts, std::size_t du, std::size_t dv, Parametrization method) {
        auto [Q, nu, nv] = grid(pts);
        return interpolate_surface(Q, nu, nv, du, dv, method);
      },
      py::arg("points"), py::arg("degree_u") = 3, py::arg("degree_v") = 3,
      py::arg("method") = Parametrization::Centripetal, "Interpolate a point grid (nu, nv, 3).");
  m.def(
      "approximate_surface",
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
  m.def("make_circle_3pt", &make_circle_3pt, py::arg("p0"), py::arg("p1"), py::arg("p2"),
        "Full circle through three points (starts at p0, runs through p1 and p2).");
  m.def("make_arc_3pt", &make_arc_3pt, py::arg("p0"), py::arg("p1"), py::arg("p2"),
        "Circular arc from p0 through p1 to p2.");
  m.def("make_arc_tangent", &make_arc_tangent, py::arg("start"), py::arg("tangent"), py::arg("end"),
        "Exact circular arc from `start`, leaving it in the direction `tangent`, to `end` (sweep "
        "between 0 and 360 degrees). Raises GeometryError if `end` lies on the tangent line.");
  m.def("make_polygon", &make_polygon, py::arg("center"), py::arg("normal"), py::arg("radius"),
        py::arg("n"), py::arg("x_axis") = Vec3{0, 0, 0},
        "Closed regular polygon with n vertices on the circle of the given radius (first vertex "
        "along x_axis).");
  m.def("make_rectangle_curve", &make_rectangle_curve, py::arg("origin"), py::arg("edge_u"),
        py::arg("edge_v"), "Closed rectangle outline (degree 1).");
  m.def("make_helix", &make_helix, py::call_guard<py::gil_scoped_release>(), py::arg("axis_point"),
        py::arg("axis"), py::arg("start"), py::arg("pitch"), py::arg("turns"),
        py::arg("end_radius") = -1.0, py::arg("left_handed") = false, py::arg("tolerance") = 1e-7,
        "Helix / conical helix / planar spiral starting at `start`, parameter = turns; C2 quintic "
        "within "
        "`tolerance` of the exact curve. end_radius < 0: cylindrical.");
  m.def("make_bilinear", &make_bilinear, py::arg("p00"), py::arg("p10"), py::arg("p01"),
        py::arg("p11"));
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
        py::call_guard<py::gil_scoped_release>(), py::arg("curve"), py::arg("end"),
        py::arg("length"), py::arg("type") = Extension::Natural,
        "Extend a curve at one end by an arc length; original parameters are kept.");
  m.def("extend", py::overload_cast<const NurbsSurface&, SurfaceSide, double, Extension>(&extend),
        py::call_guard<py::gil_scoped_release>(), py::arg("surface"), py::arg("side"),
        py::arg("length"), py::arg("type") = Extension::Natural,
        "Extend a surface across one side; original parameters are kept.");

  // ------------------------------------------------------------------ loft
  m.def("loft_parameters", &loft_parameters, py::arg("sections"),
        py::arg("method") = Parametrization::Centripetal);
  m.def("loft",
        py::overload_cast<const std::vector<NurbsCurve>&, std::size_t, Parametrization>(&loft),
        py::call_guard<py::gil_scoped_release>(), py::arg("sections"), py::arg("degree_v") = 3,
        py::arg("method") = Parametrization::Centripetal,
        "Skinned surface through section curves (u along the sections, v across them).");
  m.def("loft",
        py::overload_cast<const std::vector<NurbsCurve>&, std::size_t, const std::vector<double>&>(
            &loft),
        py::call_guard<py::gil_scoped_release>(), py::arg("sections"), py::arg("degree_v"),
        py::arg("v_params"));
  m.def(
      "make_compatible",
      [](std::vector<NurbsCurve> curves) {
        make_compatible(curves);
        return curves;
      },
      py::arg("curves"), "Common degree, knot vector and domain [0, 1] for all curves.");

  // ------------------------------------------------------------------ sweep
  py::enum_<SweepFrame>(m, "SweepFrame")
      .value("RotationMinimizing", SweepFrame::RotationMinimizing)
      .value("Fixed", SweepFrame::Fixed);
  py::class_<SweepOptions>(m, "SweepOptions")
      .def(py::pickle(
          [](const SweepOptions& o) {
            return py::make_tuple(o.frame, o.twist, o.tolerance, o.max_sections);
          },
          [](const py::tuple& t) {
            SweepOptions o;
            o.frame = t[0].cast<SweepFrame>();
            o.twist = t[1].cast<double>();
            o.tolerance = t[2].cast<double>();
            o.max_sections = t[3].cast<std::size_t>();
            return o;
          }))
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
  m.def("sweep",
        py::overload_cast<const NurbsCurve&, const NurbsCurve&, const SweepOptions&>(&sweep),
        py::call_guard<py::gil_scoped_release>(), py::arg("profile"), py::arg("trajectory"),
        py::arg("options") = SweepOptions{},
        "Sweep a profile (given at the trajectory start) along a trajectory.");
  m.def("sweep",
        py::overload_cast<const std::vector<NurbsCurve>&, const NurbsCurve&, std::vector<double>,
                          const SweepOptions&>(&sweep),
        py::call_guard<py::gil_scoped_release>(), py::arg("profiles"), py::arg("trajectory"),
        py::arg("params") = std::vector<double>{}, py::arg("options") = SweepOptions{},
        "Sweep with several profiles, each given where it crosses the trajectory (at params[k]; "
        "empty: nearest trajectory point); the shape is blended between them (loft) and placed "
        "along the trajectory.");
  m.def("rail_revolve", &rail_revolve, py::call_guard<py::gil_scoped_release>(), py::arg("profile"),
        py::arg("rail"), py::arg("axis_point"), py::arg("axis_dir"),
        py::arg("scale_height") = false, py::arg("options") = SweepOptions{},
        "Revolve a profile about an axis while its end nearer to the rail follows the rail "
        "(radial scaling about the axis; scale_height: also axial scaling about the axis "
        "point). E.g. a reflector with a non-round rim.");
  m.def("sweep2", &sweep2, py::call_guard<py::gil_scoped_release>(), py::arg("profile"),
        py::arg("rail1"), py::arg("rail2"), py::arg("options") = SweepOptions{},
        "Two-rail sweep with uniform scaling.");

  // ------------------------------------------------------------------ offsets
  m.def("offset_curve", &offset_curve, py::call_guard<py::gil_scoped_release>(), py::arg("curve"),
        py::arg("distance"), py::arg("plane_normal") = Vec3{0, 0, 1}, py::arg("tolerance") = 1e-7,
        "Planar offset c(t) + d * normalize(plane_normal x c'(t)) as a cubic B-spline fitted to "
        "the tolerance (parameters kept). Raises GeometryError at cusps (``|d|`` >= radius of "
        "curvature on the concave side), corners or non-planar curves.");
  m.def("offset_surface", &offset_surface, py::call_guard<py::gil_scoped_release>(),
        py::arg("surface"), py::arg("distance"), py::arg("tolerance") = 1e-7,
        "Offset surface S(u,v) + d * N(u,v) as a quintic tensor-product spline fitted to the "
        "tolerance (parameters kept). Raises GeometryError if d reaches a radius of curvature.");

  // ------------------------------------------------------------------ fillets
  py::class_<CurveFillet>(m, "CurveFillet", "Circular fillet between two planar curves.")
      .def_readonly("arc", &CurveFillet::arc, "exact rational arc from point_a to point_b")
      .def_readonly("center", &CurveFillet::center)
      .def_readonly("param_a", &CurveFillet::param_a, "parameter of the tangent point on a")
      .def_readonly("param_b", &CurveFillet::param_b, "parameter of the tangent point on b")
      .def_readonly("point_a", &CurveFillet::point_a)
      .def_readonly("point_b", &CurveFillet::point_b)
      .def_readonly("trimmed_a", &CurveFillet::trimmed_a,
                    "kept part of a ending at point_a (None if trim=False)")
      .def_readonly("trimmed_b", &CurveFillet::trimmed_b,
                    "kept part of b starting at point_b (None if trim=False)")
      .def("__repr__", [](const CurveFillet& f) {
        return "<CurveFillet center=(" + std::to_string(f.center.x) + ", " +
               std::to_string(f.center.y) + ", " + std::to_string(f.center.z) + ")>";
      });
  m.def("fillet_corners", &fillet_corners, py::call_guard<py::gil_scoped_release>(),
        py::arg("curve"), py::arg("radius"), py::arg("tolerance") = 1e-7,
        "Round every corner of a planar curve (polyline, polygon, ...) with exact arcs; one G1 "
        "curve.");
  m.def(
      "fillet_curves",
      [](const NurbsCurve& a, const NurbsCurve& b, double radius, const Vec3& plane_normal,
         int side_a, int side_b, bool trim, std::optional<Vec3> corner, double tolerance) {
        return fillet_curves(a, b, radius, plane_normal,
                             CurveFilletOptions{side_a, side_b, trim, corner, tolerance});
      },
      py::arg("a"), py::arg("b"), py::arg("radius"), py::arg("plane_normal") = Vec3{0, 0, 1},
      py::arg("side_a") = 0, py::arg("side_b") = 0, py::arg("trim") = true,
      py::arg("corner") = py::none(), py::arg("tolerance") = 1e-7,
      "Circular fillet of the given radius between two coplanar intersecting curves. side_a / "
      "side_b: +1 = centre to the left of the curve's direction seen from plane_normal, -1 = "
      "right, 0 = automatic (the corner between the kept parts, i.e. the longer sides of the "
      "intersection). Returns a CurveFillet (exact arc, tangent parameters, trimmed curves).");
  py::class_<CurveChamfer>(m, "CurveChamfer", "Chamfer (straight cut) between two curves.")
      .def_readonly("line", &CurveChamfer::line, "straight segment from point_a to point_b")
      .def_readonly("param_a", &CurveChamfer::param_a, "parameter of the chamfer point on a")
      .def_readonly("param_b", &CurveChamfer::param_b, "parameter of the chamfer point on b")
      .def_readonly("point_a", &CurveChamfer::point_a, "chamfer point on a (line start)")
      .def_readonly("point_b", &CurveChamfer::point_b, "chamfer point on b (line end)")
      .def_readonly("trimmed_a", &CurveChamfer::trimmed_a,
                    "kept part of a ending at point_a (None if trim=False)")
      .def_readonly("trimmed_b", &CurveChamfer::trimmed_b,
                    "kept part of b starting at point_b (None if trim=False)")
      .def("__repr__", [](const CurveChamfer& c) {
        return "<CurveChamfer length=" + std::to_string(distance(c.point_a, c.point_b)) + ">";
      });
  m.def(
      "chamfer_curves",
      [](const NurbsCurve& a, const NurbsCurve& b, double distance_a,
         std::optional<double> distance_b, bool trim, std::optional<Vec3> corner) {
        return chamfer_curves(a, b, distance_a, distance_b.value_or(distance_a),
                              CurveChamferOptions{trim, corner});
      },
      py::call_guard<py::gil_scoped_release>(), py::arg("a"), py::arg("b"), py::arg("distance_a"),
      py::arg("distance_b") = py::none(), py::arg("trim") = true, py::arg("corner") = py::none(),
      "Chamfer between two intersecting curves: a straight cut from the point of a at distance "
      "`distance_a` from the corner to the point of b at `distance_b` (default: distance_a), on "
      "the kept parts (the longer sides of the intersection, as in fillet_curves). The curves "
      "need not be planar. Returns a CurveChamfer (line, parameters, trimmed curves).");
  m.def("chamfer_corners", &chamfer_corners, py::call_guard<py::gil_scoped_release>(),
        py::arg("curve"), py::arg("distance"), py::arg("tolerance") = 1e-7,
        "Chamfer every corner of a curve (polyline, polygon, ...) with straight cuts at "
        "`distance` from the corner on both sides; one curve.");
  py::class_<SurfaceFillet>(m, "SurfaceFillet", "Rolling-ball fillet between two surfaces.")
      .def_readonly("surface", &SurfaceFillet::surface,
                    "u: arc from contact1 (u=0) to contact2 (u=1), v: spine parameter")
      .def_readonly("spine", &SurfaceFillet::spine, "path of the ball centre")
      .def_readonly("contact1", &SurfaceFillet::contact1)
      .def_readonly("contact2", &SurfaceFillet::contact2)
      .def_readonly("pcurve1", &SurfaceFillet::pcurve1, "contact1 as (u, v, 0) on s1")
      .def_readonly("pcurve2", &SurfaceFillet::pcurve2, "contact2 as (u, v, 0) on s2")
      .def("__repr__", [](const SurfaceFillet& f) {
        return "<SurfaceFillet sections=" + std::to_string(f.surface.num_v()) + ">";
      });
  m.def(
      "fillet_surfaces",
      [](const NurbsSurface& s1, const NurbsSurface& s2, double radius, int side1, int side2,
         double tolerance, double offset_tolerance, std::size_t max_sections,
         std::optional<Vec3> near) {
        py::gil_scoped_release release;
        return fillet_surfaces(
            s1, s2, radius,
            SurfaceFilletOptions{side1, side2, tolerance, offset_tolerance, max_sections, near});
      },
      py::arg("s1"), py::arg("s2"), py::arg("radius"), py::arg("side1") = 1, py::arg("side2") = 1,
      py::arg("tolerance") = 1e-6, py::arg("offset_tolerance") = 0.0,
      py::arg("max_sections") = 1025, py::arg("near") = py::none(),
      "Rolling-ball fillet of radius r between two surfaces. side1 / side2 = +1 put the ball on "
      "the side of the normal S_u x S_v of s1 / s2, -1 on the opposite side; `near` selects the "
      "spine nearest to a point if there are several. Returns a SurfaceFillet (surface, spine, "
      "contact curves, p-curves on s1 / s2).");

  // ------------------------------------------------------------------ Coons / Gordon
  m.def("coons", &coons, py::call_guard<py::gil_scoped_release>(), py::arg("c0"), py::arg("c1"),
        py::arg("d0"), py::arg("d1"), py::arg("tolerance") = 1e-7,
        "Coons patch: S(u,0)=c0, S(u,1)=c1, S(0,v)=d0, S(1,v)=d1 (orientations fixed "
        "automatically). Exact also for rational boundaries (arcs) whose corner weights agree.");
  // ------------------------------------------------------------------ surface tools
  py::class_<RebuiltCurve>(m, "RebuiltCurve")
      .def_readonly("curve", &RebuiltCurve::curve)
      .def_readonly("max_error", &RebuiltCurve::max_error,
                    "max distance to the input at equal parameters (dense samples)")
      .def("__repr__", [](const RebuiltCurve& r) {
        return "<RebuiltCurve max_error=" + std::to_string(r.max_error) + ">";
      });
  py::class_<RebuiltSurface>(m, "RebuiltSurface")
      .def_readonly("surface", &RebuiltSurface::surface)
      .def_readonly("max_error", &RebuiltSurface::max_error,
                    "max distance to the input at equal parameters (dense grid)")
      .def("__repr__", [](const RebuiltSurface& r) {
        return "<RebuiltSurface max_error=" + std::to_string(r.max_error) + ">";
      });
  m.def("rebuild_curve", &rebuild_curve, py::call_guard<py::gil_scoped_release>(), py::arg("curve"),
        py::arg("num_ctrl"), py::arg("degree") = 3,
        "Rebuild a curve as a non-rational B-spline with num_ctrl control points and uniform "
        "knots on the same domain (least squares, parametrisation kept); returns RebuiltCurve "
        "with the max error.");
  m.def("rebuild_surface", &rebuild_surface, py::call_guard<py::gil_scoped_release>(),
        py::arg("surface"), py::arg("num_u"), py::arg("num_v"), py::arg("degree_u") = 3,
        py::arg("degree_v") = 3,
        "Rebuild a surface with num_u x num_v control points and uniform knots on the same "
        "domain; returns RebuiltSurface with the max error.");
  m.def("extrude_tapered", &extrude_tapered, py::call_guard<py::gil_scoped_release>(),
        py::arg("profile"), py::arg("direction"), py::arg("angle"), py::arg("tolerance") = 1e-7,
        "Tapered (draft) extrusion of a planar profile normal to `direction` (its length is the "
        "height): the top is the profile offset by height*tan(angle), outwards for angle > 0. "
        "Corners are mitred; returns one ruled surface per smooth piece of the profile.");
  m.def("ribbon", &ribbon, py::call_guard<py::gil_scoped_release>(), py::arg("curve"),
        py::arg("width"), py::arg("plane_normal") = Vec3{0, 0, 1}, py::arg("tolerance") = 1e-7,
        "Planar ribbon: ruled surfaces between a curve and its offset by `width` (left of the "
        "curve seen from the normal's tip), mitred at corners; one surface per smooth piece.");
  m.def("extrude_to_point", &extrude_to_point, py::arg("curve"), py::arg("apex"),
        "Cone over a curve: ruled surface from the curve (v = 0) to the apex (v = 1); exact, "
        "also for rational curves.");
  m.def("triangle_patch", &triangle_patch, py::call_guard<py::gil_scoped_release>(), py::arg("c0"),
        py::arg("c1"), py::arg("c2"), py::arg("tolerance") = 1e-7,
        "Three-sided patch bounded by three curves forming a loop: a Coons patch with c0 as one "
        "side and the opposite corner as a degenerate side (v = 1).");
  m.def("gordon", &gordon, py::call_guard<py::gil_scoped_release>(), py::arg("u_curves"),
        py::arg("v_curves"), py::arg("tolerance") = 1e-7,
        "Gordon surface through a consistent curve network.");

  // ------------------------------------------------------------------ optical surfaces
  py::class_<Asphere>(
      m, "Asphere",
      "Rotationally symmetric sag: conic + even (r^4, r^6, ...) + odd (r, r^2, ...) terms.")
      .def(py::pickle(
          [](const Asphere& a) { return py::make_tuple(a.curvature, a.conic, a.even, a.odd); },
          [](const py::tuple& t) {
            Asphere a;
            a.curvature = t[0].cast<double>();
            a.conic = t[1].cast<double>();
            a.even = t[2].cast<std::vector<double>>();
            a.odd = t[3].cast<std::vector<double>>();
            return a;
          }))
      .def(py::init([](double curvature, double conic, std::vector<double> even,
                       std::vector<double> odd) {
             return Asphere{curvature, conic, std::move(even), std::move(odd)};
           }),
           py::arg("curvature") = 0.0, py::arg("conic") = 0.0,
           py::arg("even") = std::vector<double>{}, py::arg("odd") = std::vector<double>{})
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
  m.def(
      "zernike_nm",
      [](std::size_t j, ZernikeIndexing ix) {
        int n = 0, mm = 0;
        zernike_nm(j, ix, n, mm);
        return py::make_tuple(n, mm);
      },
      py::arg("j"), py::arg("indexing"));
  py::class_<ZernikeSag>(m, "ZernikeSag")
      .def(py::pickle(
          [](const ZernikeSag& z) {
            return py::make_tuple(z.base, z.norm_radius, z.indexing, z.coeffs);
          },
          [](const py::tuple& t) {
            ZernikeSag z;
            z.base = t[0].cast<Asphere>();
            z.norm_radius = t[1].cast<double>();
            z.indexing = t[2].cast<ZernikeIndexing>();
            z.coeffs = t[3].cast<std::vector<double>>();
            return z;
          }))
      .def(
          py::init([](Asphere base, double norm_radius, ZernikeIndexing ix, std::vector<double> c) {
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
      .def(py::pickle(
          [](const XYPolynomialSag& x) {
            py::list terms;
            for (const auto& term : x.terms) terms.append(py::make_tuple(term.i, term.j, term.c));
            return py::make_tuple(x.base, x.norm_radius, terms);
          },
          [](const py::tuple& t) {
            XYPolynomialSag x;
            x.base = t[0].cast<Asphere>();
            x.norm_radius = t[1].cast<double>();
            for (const auto& item : t[2].cast<py::list>()) {
              const auto e = item.cast<py::tuple>();
              x.terms.push_back({e[0].cast<int>(), e[1].cast<int>(), e[2].cast<double>()});
            }
            return x;
          }))
      .def(py::init([](Asphere base, double norm_radius,
                       const std::vector<std::tuple<int, int, double>>& terms) {
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
        py::arg("y1"), py::arg("degree") = 3, py::arg("tolerance") = 1e-7,
        py::arg("max_spans") = 256,
        "Graph surface z = f(x, y) over a rectangle, fitted to the tolerance.");
  m.def("zernike_surface", &zernike_surface, py::call_guard<py::gil_scoped_release>(),
        py::arg("zernike"), py::arg("half_width") = 0.0, py::arg("tolerance") = 1e-7);
  m.def("xy_polynomial_surface", &xy_polynomial_surface, py::call_guard<py::gil_scoped_release>(),
        py::arg("poly"), py::arg("half_width") = 0.0, py::arg("tolerance") = 1e-7);

  // ------------------------------------------------------------------ trimmed surfaces
  py::class_<TrimLoop>(m, "TrimLoop",
                       "Closed loop of parameter-space curves (x = u, y = v, z = 0).")
      .def(py::pickle(
          [](const TrimLoop& l) { return py::make_tuple(l.curves); },
          [](const py::tuple& t) { return TrimLoop{t[0].cast<std::vector<NurbsCurve>>()}; }))
      .def(py::init([](std::vector<NurbsCurve> curves) { return TrimLoop{std::move(curves)}; }),
           py::arg("curves"))
      .def_readonly("curves", &TrimLoop::curves)
      .def("signed_area", &TrimLoop::signed_area)
      .def("reversed", &TrimLoop::reversed)
      .def("polyline", &TrimLoop::polyline, py::arg("samples_per_span") = 16);
  m.def("rectangle_loop", &rectangle_loop, py::arg("u0"), py::arg("u1"), py::arg("v0"),
        py::arg("v1"));
  m.def("circle_loop", &circle_loop, py::arg("cu"), py::arg("cv"), py::arg("r"));
  m.def("polygon_loop", &polygon_loop, py::arg("points"));
  py::class_<TrimmedSurface>(m, "TrimmedSurface")
      .def(py::pickle(
          [](const TrimmedSurface& f) {
            return py::make_tuple(f.surface(), f.outer(), f.holes(), f.label());
          },
          [](const py::tuple& t) {
            return TrimmedSurface(t[0].cast<NurbsSurface>(), t[1].cast<TrimLoop>(),
                                  t[2].cast<std::vector<TrimLoop>>())
                .with_label(t[3].cast<std::string>());
          }))
      .def(py::init<NurbsSurface>(), py::arg("surface"))
      .def(py::init<NurbsSurface, TrimLoop, std::vector<TrimLoop>, const Tolerance&>(),
           py::arg("surface"), py::arg("outer"), py::arg("holes") = std::vector<TrimLoop>{},
           py::arg("tolerance") = Tolerance{})
      .def_property_readonly("surface", &TrimmedSurface::surface)
      .def_property_readonly("outer", &TrimmedSurface::outer)
      .def_property_readonly("holes", &TrimmedSurface::holes)
      .def_property_readonly("is_untrimmed", &TrimmedSurface::is_untrimmed)
      .def("contains", &TrimmedSurface::contains, py::arg("u"), py::arg("v"))
      .def("area", &TrimmedSurface::area, py::call_guard<py::gil_scoped_release>(),
           py::arg("rel_tol") = 1e-9)
      .def("boundary_curves", &TrimmedSurface::boundary_curves,
           py::call_guard<py::gil_scoped_release>(), py::arg("tolerance") = 1e-7)
      .def("with_hole", &TrimmedSurface::with_hole, py::arg("hole"))
      .def("flipped", &TrimmedSurface::flipped,
           "The same face with opposite orientation (normal flipped).")
      .def("transformed", &TrimmedSurface::transformed, py::arg("transform"),
           "The face mapped by an affine transform (loops and label kept; reflections keep the "
           "normal by flipping the face).")
      .def("__rmul__", [](const TrimmedSurface& f, const Transform& t) { return f.transformed(t); })
      .def_property_readonly(
          "label", &TrimmedSurface::label,
          "Face label (boundary condition / physical group); kept by Booleans, transforms, STEP.")
      .def("with_label", &TrimmedSurface::with_label, py::arg("label"),
           "The same face with the given label.")
      .def("__repr__", [](const TrimmedSurface& t) {
        return "<TrimmedSurface holes=" + std::to_string(t.holes().size()) + ">";
      });
  m.def("split",
        py::overload_cast<const NurbsSurface&, const NurbsSurface&, const Tolerance&>(&split),
        py::call_guard<py::gil_scoped_release>(), py::arg("surface"), py::arg("cutter"),
        py::arg("tolerance") = Tolerance{},
        "Split a surface along its intersection curves with another surface into trimmed regions.");
  m.def(
      "split",
      py::overload_cast<const NurbsSurface&, const std::vector<NurbsCurve>&, const Tolerance&>(
          &split),
      py::call_guard<py::gil_scoped_release>(), py::arg("surface"), py::arg("pcurves"),
      py::arg("tolerance") = Tolerance{},
      "Split a surface into trimmed regions along parameter-space curves (u, v, 0). Closed curves "
      "become holes and new regions; open curves run from boundary to boundary and may cross.");
  m.def(
      "split",
      py::overload_cast<const TrimmedSurface&, const std::vector<NurbsCurve>&, const Tolerance&>(
          &split),
      py::call_guard<py::gil_scoped_release>(), py::arg("face"), py::arg("cuts"),
      py::arg("tolerance") = Tolerance{},
      "Split a trimmed face along parameter-space curves (general planar arrangement): cuts are "
      "clipped to the face, may end on the trim loops or on other cuts and may cross each other.");
  m.def("split", py::overload_cast<const Shell&, const Plane&, double, const std::string&>(&split),
        py::call_guard<py::gil_scoped_release>(), py::arg("shell"), py::arg("plane"),
        py::arg("tolerance") = 1e-6, py::arg("cut_label") = "cut",
        "Split a closed shell by a plane -> (below, above): two closed shells, below opposite to "
        "the plane normal. The new faces in the plane get `cut_label` (e.g. for a symmetry "
        "boundary condition); a part is an empty shell if the plane misses the solid.");
  m.def("cap", &cap, py::call_guard<py::gil_scoped_release>(), py::arg("shell"),
        py::arg("tolerance") = 1e-6, py::arg("label") = "cap",
        "Close the planar openings of a shell with planar faces trimmed exactly by the opening's "
        "edges (coplanar nested openings: faces with holes). The new faces get `label`; the "
        "result is sewn and oriented outwards. Raises GeometryError for non-planar openings.");
  m.def("thicken", &thicken, py::call_guard<py::gil_scoped_release>(), py::arg("shell"),
        py::arg("thickness"), py::arg("symmetric") = false, py::arg("tolerance") = 1e-6,
        py::arg("wall_label") = "wall",
        "Give an open shell a wall thickness: faces offset along S_u x S_v (negative: against it; "
        "symmetric: +-thickness/2), free edges closed by ruled side walls labelled `wall_label`, "
        "offset faces labelled '<label>_offset'. The faces must join smoothly; raises "
        "GeometryError at creases or where the thickness reaches a radius of curvature.");
  m.def("make_pipe", py::overload_cast<const NurbsCurve&, double, double, double>(&make_pipe),
        py::call_guard<py::gil_scoped_release>(), py::arg("path"), py::arg("radius"),
        py::arg("end_radius") = -1.0, py::arg("tolerance") = 1e-6,
        "Solid pipe along a path: circular cross-section (radius varying linearly to "
        "`end_radius` if given), planar end faces. Faces labelled 'wall', 'inlet', 'outlet'; a "
        "closed path gives a ring.");
  m.def("make_pipe",
        py::overload_cast<const NurbsCurve&, const std::function<double(double)>&, double>(
            &make_pipe),
        py::arg("path"), py::arg("radius"), py::arg("tolerance") = 1e-6,
        "Solid pipe whose radius is the function radius(t) of the normalised path parameter t in "
        "[0, 1] (smooth and positive).");
  m.def("pipe_surface", &pipe_surface, py::call_guard<py::gil_scoped_release>(), py::arg("path"),
        py::arg("radius"), py::arg("options") = SweepOptions{},
        "Pipe wall: a circle of radius radius(t) (t in [0, 1] along the path) swept with the "
        "rotation-minimising frame; u around the circle, v along the path.");
  m.def("remove_holes", &remove_holes, py::call_guard<py::gil_scoped_release>(), py::arg("shell"),
        py::arg("max_size"), py::arg("tolerance") = 1e-6,
        "Defeaturing: fill the holes (through or blind) of a closed shell whose loops are smaller "
        "than max_size across; the faces of the removed features are deleted.");
  m.def("edges_between", &edges_between, py::arg("shell"), py::arg("label_a"), py::arg("label_b"),
        "Indices of the edges shared by a face labelled `label_a` and one labelled `label_b`.");
  const auto blend_options = [](const std::string& label, double tolerance,
                                const std::vector<std::pair<Vec3, double>>& sizes) {
    EdgeBlendOptions o;
    o.label = label;
    o.tolerance = tolerance;
    for (const auto& [p, v] : sizes) o.sizes.push_back({p, v});
    return o;
  };
  m.def(
      "fillet_edges",
      [blend_options](const Shell& shell, const std::vector<std::size_t>& edges, double radius,
                      const std::string& label, double tolerance,
                      const std::vector<std::pair<Vec3, double>>& sizes) {
        return fillet_edges(shell, edges, radius, blend_options(label, tolerance, sizes));
      },
      py::call_guard<py::gil_scoped_release>(), py::arg("shell"), py::arg("edges"),
      py::arg("radius"), py::arg("label") = "fillet", py::arg("tolerance") = 0.0,
      py::arg("sizes") = std::vector<std::pair<Vec3, double>>{},
      "Round edges of a closed shell with a rolling-ball fillet (convex edges lose, concave "
      "edges gain material). Edges between the same two faces are filleted as one chain. Where "
      "chains meet at a vertex of three faces: tangent chains continue each other, two chains "
      "are mitred (trimmed by each other), three chains get a spherical corner patch; a chain "
      "ending alone is trimmed by the third face. Variable radius: `sizes` = [(point, radius), "
      "...] at points on the selected edges, smooth (C2) in between, constant beyond; chains "
      "without points use `radius`. New faces are labelled `label`.");
  m.def(
      "chamfer_edges",
      [blend_options](const Shell& shell, const std::vector<std::size_t>& edges, double distance,
                      const std::string& label, double tolerance,
                      const std::vector<std::pair<Vec3, double>>& sizes) {
        return chamfer_edges(shell, edges, distance, blend_options(label, tolerance, sizes));
      },
      py::call_guard<py::gil_scoped_release>(), py::arg("shell"), py::arg("edges"),
      py::arg("distance"), py::arg("label") = "chamfer", py::arg("tolerance") = 0.0,
      py::arg("sizes") = std::vector<std::pair<Vec3, double>>{},
      "Chamfer edges of a closed shell: a ruled cut between the curves at (geodesic) `distance` "
      "from the edge on both faces. Grouping, vertices, variable `sizes` (point, distance) and "
      "restrictions as for fillet_edges; a corner of three chamfered edges (between planar "
      "faces) gets a planar triangle.");
  m.def("merge_faces", &merge_faces, py::call_guard<py::gil_scoped_release>(), py::arg("shell"),
        py::arg("tolerance") = 1e-6,
        "Defeaturing: merge adjacent faces with the same label lying on the same surface (pieces "
        "of one split face, or coplanar planar faces with the same normal) into one face.");
  m.def("hollow", &hollow, py::call_guard<py::gil_scoped_release>(), py::arg("shell"),
        py::arg("thickness"), py::arg("open_labels") = std::vector<std::string>{},
        py::arg("tolerance") = 1e-6,
        "Hollow a smooth closed shell to a wall of the given thickness: an inner void, or - with "
        "open_labels - removing the faces with these labels and thickening the rest inwards.");
  m.def("blend_curves", &blend_curves, py::call_guard<py::gil_scoped_release>(), py::arg("a"),
        py::arg("end0"), py::arg("b"), py::arg("end1"), py::arg("continuity") = 2,
        py::arg("tension") = 1.0,
        "Smooth transition from an end of curve a to an end of curve b: G1 (cubic) or G2 "
        "(quintic, same curvature vectors). The blend leaves a and enters b; |C'| at the ends is "
        "tension x chord.");
  m.def("blend_surfaces", &blend_surfaces, py::call_guard<py::gil_scoped_release>(), py::arg("s0"),
        py::arg("side0"), py::arg("s1"), py::arg("side1"), py::arg("continuity") = 2,
        py::arg("tension") = 1.0, py::arg("tolerance") = 1e-7,
        "Blend surface between a side of s0 and a side of s1 continuing both with G0 (ruled), "
        "G1 (cubic across) or G2 (quintic across, same normal curvature). u runs along the edge "
        "of s0, v from s0 (v = 0) to s1 (v = 1); the edge of s1 is mapped linearly onto it in "
        "the orientation pairing the nearer points. The cross speed is tension x the mean edge "
        "distance. Exact when the weights do not vary across the sides, otherwise the boundary "
        "data are interpolated to tolerance.");
  m.def("match_surface", &match_surface, py::call_guard<py::gil_scoped_release>(),
        py::arg("surface"), py::arg("side"), py::arg("target"), py::arg("target_side"),
        py::arg("continuity") = 2, py::arg("tolerance") = 1e-7,
        "Modify a surface near one side so that it joins the target's side with G0, G1 or G2 "
        "continuity: the first continuity + 1 control rows are replaced (after refinement), the "
        "opposite side and its cross derivatives are kept. Same parameter directions and "
        "domain.");
  py::class_<ProjectedCurve>(m, "ProjectedCurve", "A curve on a surface: 3D curve and p-curve.")
      .def_readonly("curve", &ProjectedCurve::curve, "3D curve on the surface.")
      .def_readonly("pcurve", &ProjectedCurve::pcurve,
                    "(u, v, 0) in the surface's parameter domain (same parameter as curve).")
      .def("__iter__",
           [](const ProjectedCurve& p) { return py::iter(py::make_tuple(p.curve, p.pcurve)); });
  m.def("project_curve", &project_curve, py::call_guard<py::gil_scoped_release>(), py::arg("curve"),
        py::arg("surface"), py::arg("direction") = Vec3{0, 0, 0}, py::arg("tolerance") = 1e-6,
        "Project a curve onto a surface along its normal (direction zero) or along `direction`; "
        "returns ProjectedCurve(curve, pcurve) within `tolerance` of the exact projection.");
  m.def(
      "project_to_plane",
      [](const NurbsCurve& c, const Plane& plane, const Vec3& direction) {
        return project_to_plane(c, plane, direction);
      },
      py::arg("curve"), py::arg("plane"), py::arg("direction") = Vec3{0, 0, 0},
      "Exact parallel projection of a curve onto a plane along `direction` (zero: along the "
      "plane normal, i.e. orthogonal projection).");
  m.def(
      "project_to_plane",
      [](const NurbsSurface& s, const Plane& plane, const Vec3& direction) {
        return project_to_plane(s, plane, direction);
      },
      py::arg("surface"), py::arg("plane"), py::arg("direction") = Vec3{0, 0, 0},
      "Exact parallel projection of a surface onto a plane.");
  m.def(
      "interpolate_on_surface",
      [](const NurbsSurface& s, const Arr& pts, std::size_t degree, bool uv, double tolerance) {
        auto P = from_numpy(pts);
        py::gil_scoped_release release;
        return interpolate_on_surface(s, P, degree, uv, tolerance);
      },
      py::arg("surface"), py::arg("points"), py::arg("degree") = 3, py::arg("uv") = false,
      py::arg("tolerance") = 1e-6,
      "Curve on a surface through points (n, 3): projected onto the surface (or (u, v, 0) "
      "parameter points if uv=True), interpolated in the parameter domain; returns "
      "ProjectedCurve(curve, pcurve) with the 3D curve within `tolerance` of S(pcurve(t)).");
  m.def(
      "offset_on_surface",
      [](const NurbsSurface& s, const py::object& c, double distance, double tolerance) {
        const NurbsCurve pc = py::isinstance<ProjectedCurve>(c)
                                  ? c.cast<const ProjectedCurve&>().pcurve
                                  : c.cast<NurbsCurve>();
        py::gil_scoped_release release;
        return offset_on_surface(s, pc, distance, tolerance);
      },
      py::arg("surface"), py::arg("curve"), py::arg("distance"), py::arg("tolerance") = 1e-6,
      "Geodesic offset of a curve within a surface: `curve` is its p-curve (u, v, 0) or a "
      "ProjectedCurve; positive distances go to the left of the curve seen from the normal "
      "S_u x S_v. Returns ProjectedCurve(curve, pcurve) within `tolerance`.");
  m.def("tween_curves", &tween_curves, py::call_guard<py::gil_scoped_release>(), py::arg("a"),
        py::arg("b"), py::arg("count"),
        "`count` intermediate curves between a and b (made compatible, homogeneous control "
        "points interpolated linearly at s = k / (count + 1)). Both should run in the same "
        "direction.");
  py::class_<TangentLine>(m, "TangentLine", "Line segment tangent to curves.")
      .def_readonly("line", &TangentLine::line, "segment from point_a to point_b")
      .def_readonly("point_a", &TangentLine::point_a,
                    "line start: tangent point on the first curve, or the given point")
      .def_readonly("point_b", &TangentLine::point_b, "line end: tangent point on the second curve")
      .def_readonly("param_a", &TangentLine::param_a,
                    "tangent parameter on the first curve (nan if it is a point)")
      .def_readonly("param_b", &TangentLine::param_b, "tangent parameter on the second curve")
      .def("__repr__", [](const TangentLine& l) {
        return "<TangentLine length=" + std::to_string(distance(l.point_a, l.point_b)) + ">";
      });
  m.def(
      "tangent_lines",
      [](const NurbsCurve& a, const NurbsCurve& b, const Vec3& n) {
        return tangent_lines(a, b, n);
      },
      py::call_guard<py::gil_scoped_release>(), py::arg("a"), py::arg("b"),
      py::arg("plane_normal") = Vec3{0, 0, 1},
      "Common tangent lines of two coplanar curves (from a to b), or - with a point as `a` - the "
      "lines from the point tangent to curve `b`. Returns a list of TangentLine.");
  m.def(
      "tangent_lines",
      [](const Vec3& p, const NurbsCurve& b, const Vec3& n) { return tangent_lines(p, b, n); },
      py::call_guard<py::gil_scoped_release>(), py::arg("a"), py::arg("b"),
      py::arg("plane_normal") = Vec3{0, 0, 1});
  m.def(
      "region_boolean",
      [](const py::object& a, const py::object& b, BooleanOp op, double tolerance) {
        const Region ra = to_region(a), rb = to_region(b);
        py::gil_scoped_release release;
        return region_boolean(ra, rb, op, tolerance);
      },
      py::arg("a"), py::arg("b"), py::arg("op"), py::arg("tolerance") = 1e-9,
      "Union / intersection / difference of planar regions in the xy plane. A region is a "
      "PlanarSection, or a list of loops: closed curves, lists of curves running head to tail, or "
      "SectionLoops (nested loops: holes). Returns the result as a list of loops (lists of exact "
      "curves; outer counter-clockwise, holes clockwise).");
  m.def(
      "region_area", [](const py::object& r) { return region_area(to_region(r)); },
      py::arg("region"),
      "Area of a planar region (outer loops minus holes; Green's theorem on the exact curves).");
  m.def("integrate", &integrate, py::arg("f"), py::arg("a"), py::arg("b"),
        py::arg("rel_tol") = kIntegrationRelTol, py::arg("breaks") = std::vector<double>{},
        py::arg("abs_tol") = 0.0, "Globally adaptive Gauss-Kronrod integral of f over [a, b].");

  // ------------------------------------------------------------------ B-rep shells
  py::class_<Coedge>(m, "Coedge")
      .def_readonly("face", &Coedge::face)
      .def_readonly("loop", &Coedge::loop)
      .def_readonly("index", &Coedge::index)
      .def_readonly("reversed", &Coedge::reversed);
  py::class_<ShellEdge>(m, "ShellEdge")
      .def_readonly("curve", &ShellEdge::curve)
      .def_readonly("v0", &ShellEdge::v0)
      .def_readonly("v1", &ShellEdge::v1)
      .def_readonly("uses", &ShellEdge::uses);
  py::class_<MassProperties>(m, "MassProperties",
                             "Volume, centroid and inertia of a solid (unit density).")
      .def_readonly("volume", &MassProperties::volume)
      .def_readonly("centroid", &MassProperties::centroid)
      .def_property_readonly(
          "inertia",
          [](const MassProperties& mp) {
            py::array_t<double> a({py::ssize_t{3}, py::ssize_t{3}});
            auto r = a.mutable_unchecked<2>();
            for (py::ssize_t i = 0; i < 3; ++i)
              for (py::ssize_t j = 0; j < 3; ++j)
                r(i, j) = mp.inertia[static_cast<std::size_t>(3 * i + j)];
            return a;
          },
          "Inertia tensor about the centroid (3 x 3): diagonal moments, off-diagonal -products.")
      .def_property_readonly("principal_moments",
                             [](const MassProperties& mp) {
                               return py::make_tuple(mp.principal_moments[0],
                                                     mp.principal_moments[1],
                                                     mp.principal_moments[2]);
                             })
      .def_property_readonly(
          "principal_axes",
          [](const MassProperties& mp) {
            return to_numpy({mp.principal_axes[0], mp.principal_axes[1], mp.principal_axes[2]});
          },
          "Unit axes as rows, in the order of principal_moments.")
      .def("__repr__", [](const MassProperties& mp) {
        return "MassProperties(volume=" + std::to_string(mp.volume) + ", centroid=(" +
               std::to_string(mp.centroid.x) + ", " + std::to_string(mp.centroid.y) + ", " +
               std::to_string(mp.centroid.z) + "))";
      });
  py::class_<Shell>(m, "Shell", "Trimmed faces sewn along shared edges.")
      .def(py::pickle([](const Shell& s) { return py::make_tuple(s.faces(), s.tolerance()); },
                      [](const py::tuple& t) {
                        return Shell::sew(t[0].cast<std::vector<TrimmedSurface>>(),
                                          t[1].cast<double>());
                      }))
      .def_static("sew", &Shell::sew, py::call_guard<py::gil_scoped_release>(), py::arg("faces"),
                  py::arg("tolerance") = 1e-6)
      .def_property_readonly("faces", &Shell::faces)
      .def_property_readonly("vertices", &Shell::vertices)
      .def_property_readonly("edges", &Shell::edges)
      .def("free_edges", &Shell::free_edges)
      .def("non_manifold_edges", &Shell::non_manifold_edges)
      .def_property_readonly("is_closed", &Shell::is_closed)
      .def_property_readonly("is_consistently_oriented", &Shell::is_consistently_oriented)
      .def("components", &Shell::components, py::call_guard<py::gil_scoped_release>(),
           "Connected components as separate shells (outer boundary, void boundaries, lumps).")
      .def("oriented", &Shell::oriented, py::call_guard<py::gil_scoped_release>())
      .def("area", &Shell::area, py::call_guard<py::gil_scoped_release>(),
           py::arg("rel_tol") = 1e-8)
      .def("volume", &Shell::volume, py::call_guard<py::gil_scoped_release>(),
           py::arg("rel_tol") = 1e-8)
      .def("mass_properties", &Shell::mass_properties, py::call_guard<py::gil_scoped_release>(),
           py::arg("rel_tol") = 1e-8,
           "Volume, centroid and inertia tensor (about the centroid, unit density) of the solid.")
      .def("with_face_labels", &Shell::with_face_labels, py::arg("labels"),
           "The same shell with one label per face (they survive Booleans, transforms and STEP).")
      .def("faces_labeled", &Shell::faces_labeled, py::arg("label"),
           "Indices of the faces with this label.")
      .def_property_readonly("face_labels",
                             [](const Shell& s) {
                               std::vector<std::string> l;
                               for (const auto& f : s.faces()) l.push_back(f.label());
                               return l;
                             })
      .def("transformed", &Shell::transformed, py::call_guard<py::gil_scoped_release>(),
           py::arg("transform"),
           "The shell mapped by an affine transform (mirrors keep the normals outward).")
      .def("__rmul__", [](const Shell& s, const Transform& t) { return s.transformed(t); })
      .def("__repr__", [](const Shell& s) {
        return "<Shell faces=" + std::to_string(s.faces().size()) +
               " edges=" + std::to_string(s.edges().size()) +
               (s.is_closed() ? " closed>" : " open>");
      });

  // ------------------------------------------------------------------ Booleans
  py::enum_<BooleanOp>(m, "BooleanOp")
      .value("Union", BooleanOp::Union)
      .value("Intersection", BooleanOp::Intersection)
      .value("Difference", BooleanOp::Difference);
  m.def("point_in_shell", &point_in_shell, py::call_guard<py::gil_scoped_release>(),
        py::arg("shell"), py::arg("point"), py::arg("tolerance") = 1e-6,
        "True if the point lies inside the closed shell (ray casting).");
  // closed solid primitives
  m.def("make_box", &make_box, py::arg("corner"), py::arg("size"),
        "Axis-aligned box (closed Shell, outward normals) with a corner and edge lengths.");
  m.def("make_sphere", &make_sphere, py::arg("center"), py::arg("radius"));
  m.def("make_cylinder", &make_cylinder, py::arg("base"), py::arg("axis"), py::arg("radius"),
        py::arg("height"), "Cylinder with the base disc centred at `base` along `axis`.");
  m.def("make_cone", &make_cone, py::arg("base"), py::arg("axis"), py::arg("r0"), py::arg("r1"),
        py::arg("height"),
        "Frustum: radius r0 at the base, r1 at the top; one radius may be 0 (apex).");
  m.def("make_torus", &make_torus, py::arg("center"), py::arg("axis"), py::arg("major_radius"),
        py::arg("minor_radius"));
  m.def("boolean", &boolean, py::call_guard<py::gil_scoped_release>(), py::arg("a"), py::arg("b"),
        py::arg("op"), py::arg("tolerance") = 1e-6,
        "Boolean union / intersection / difference of two closed shells.");

  // ------------------------------------------------------------------ STEP export
  const auto step_unit = [](const std::string& unit) {
    if (unit == "mm") return StepLengthUnit::Millimetre;
    if (unit == "cm") return StepLengthUnit::Centimetre;
    if (unit == "m") return StepLengthUnit::Metre;
    if (unit == "inch" || unit == "in") return StepLengthUnit::Inch;
    if (unit == "um" || unit == "\u00b5m") return StepLengthUnit::Micrometre;
    throw py::value_error("unit must be 'mm', 'cm', 'm', 'inch' or 'um'");
  };
  {
    const auto shell_list = [](const py::object& shells) {
      if (py::isinstance<Shell>(shells)) return std::vector<Shell>{shells.cast<Shell>()};
      return shells.cast<std::vector<Shell>>();
    };
    const auto step_options =
        [step_unit](const std::string& unit, const std::string& name, const std::string& author,
                    const std::string& organization, const std::string& timestamp,
                    double uncertainty, const std::vector<std::string>& names) {
          StepOptions o;
          o.unit = step_unit(unit);
          o.name = name;
          o.author = author;
          o.organization = organization;
          o.timestamp = timestamp;
          o.uncertainty = uncertainty;
          o.shell_names = names;
          return o;
        };
    m.def(
        "write_step",
        [shell_list, step_options](const py::object& path, const py::object& shells,
                                   const std::string& unit, const std::string& name,
                                   const std::string& author, const std::string& organization,
                                   const std::string& timestamp, double uncertainty,
                                   const std::vector<std::string>& names) {
          const auto os = py::module_::import("os");
          write_step(os.attr("fspath")(path).cast<std::string>(), shell_list(shells),
                     step_options(unit, name, author, organization, timestamp, uncertainty, names));
        },
        py::arg("path"), py::arg("shells"), py::arg("unit") = "mm", py::arg("name") = "nurbforge",
        py::arg("author") = "", py::arg("organization") = "", py::arg("timestamp") = "",
        py::arg("uncertainty") = 1e-6, py::arg("names") = std::vector<std::string>{},
        "Write a Shell or a list of Shells to a STEP file (AP214, native writer). Closed shells "
        "become "
        "solids, open shells surface models. `unit` ('mm', 'cm', 'm', 'inch', 'um') declares what "
        "one model "
        "unit means; coordinates are written unchanged. `names`: one name per shell (solid names, "
        "read back "
        "by read_step_named).");
    m.def(
        "to_step",
        [shell_list, step_options](const py::object& shells, const std::string& unit,
                                   const std::string& name, const std::string& author,
                                   const std::string& organization, const std::string& timestamp,
                                   double uncertainty, const std::vector<std::string>& names) {
          return to_step(shell_list(shells), step_options(unit, name, author, organization,
                                                          timestamp, uncertainty, names));
        },
        py::arg("shells"), py::arg("unit") = "mm", py::arg("name") = "nurbforge",
        py::arg("author") = "", py::arg("organization") = "", py::arg("timestamp") = "",
        py::arg("uncertainty") = 1e-6, py::arg("names") = std::vector<std::string>{},
        "The STEP file contents (see write_step) as a string.");
  }

  // ------------------------------------------------------------------ STEP import
  {
    // runs the reader with the given options; skipped items are reported as Python warnings
    const auto read = [step_unit](auto&& reader, const std::string& unit, double tolerance,
                                  bool skip_unsupported) {
      StepReadOptions o;
      o.unit = step_unit(unit);
      o.tolerance = tolerance;
      o.skip_unsupported = skip_unsupported;
      std::vector<std::string> warnings;
      o.warnings = &warnings;
      decltype(reader(o)) shells;
      {
        py::gil_scoped_release release;
        shells = reader(o);
      }
      if (!warnings.empty()) {
        const auto warn = py::module_::import("warnings").attr("warn");
        for (const auto& w : warnings)
          warn(w, py::module_::import("builtins").attr("UserWarning"), 2);
      }
      return shells;
    };
    m.def(
        "read_step",
        [read](const py::object& path, const std::string& unit, double tolerance,
               bool skip_unsupported) {
          const std::string p = py::module_::import("os").attr("fspath")(path).cast<std::string>();
          return read([&](const StepReadOptions& o) { return read_step(p, o); }, unit, tolerance,
                      skip_unsupported);
        },
        py::arg("path"), py::kw_only(), py::arg("unit") = "mm", py::arg("tolerance") = 1e-6,
        py::arg("skip_unsupported") = false,
        "Read the B-rep solids and surface models of a STEP file (native reader) as a list of "
        "Shells, "
        "converted to `unit` ('mm', 'cm', 'm', 'inch', 'um'). Analytic surfaces (plane, cylinder, "
        "cone, sphere, "
        "torus) become exact NURBS patches; trim loops are rebuilt in parameter space. `tolerance` "
        "is the "
        "sewing tolerance. Unsupported entities raise GeometryError unless `skip_unsupported` (the "
        "faces "
        "are then skipped with a warning).");
    m.def(
        "from_step",
        [read](const std::string& text, const std::string& unit, double tolerance,
               bool skip_unsupported) {
          return read([&](const StepReadOptions& o) { return from_step(text, o); }, unit, tolerance,
                      skip_unsupported);
        },
        py::arg("text"), py::kw_only(), py::arg("unit") = "mm", py::arg("tolerance") = 1e-6,
        py::arg("skip_unsupported") = false,
        "Read shells from STEP file contents (see read_step).");
    const auto named = [](const std::vector<NamedShell>& v) {
      py::list out;
      for (const auto& n : v) out.append(py::make_tuple(n.name, n.shell));
      return out;
    };
    m.def(
        "read_step_named",
        [read, named](const py::object& path, const std::string& unit, double tolerance,
                      bool skip_unsupported) {
          const std::string p = py::module_::import("os").attr("fspath")(path).cast<std::string>();
          return named(read([&](const StepReadOptions& o) { return read_step_named(p, o); }, unit,
                            tolerance, skip_unsupported));
        },
        py::arg("path"), py::kw_only(), py::arg("unit") = "mm", py::arg("tolerance") = 1e-6,
        py::arg("skip_unsupported") = false,
        "As read_step, as a list of (name, Shell): the solid name (write_step(..., names=...)), "
        "else the "
        "product name.");
    m.def(
        "from_step_named",
        [read, named](const std::string& text, const std::string& unit, double tolerance,
                      bool skip_unsupported) {
          return named(read([&](const StepReadOptions& o) { return from_step_named(text, o); },
                            unit, tolerance, skip_unsupported));
        },
        py::arg("text"), py::kw_only(), py::arg("unit") = "mm", py::arg("tolerance") = 1e-6,
        py::arg("skip_unsupported") = false, "As from_step, as a list of (name, Shell).");
  }

  // ------------------------------------------------------------------ patterns
  m.def("linear_pattern", &linear_pattern, py::arg("step"), py::arg("count"));
  m.def("grid_pattern", &grid_pattern, py::arg("step_a"), py::arg("count_a"), py::arg("step_b"),
        py::arg("count_b"));
  m.def("polar_pattern", &polar_pattern, py::arg("axis_point"), py::arg("axis_dir"),
        py::arg("count"), py::arg("total_angle"));
}
