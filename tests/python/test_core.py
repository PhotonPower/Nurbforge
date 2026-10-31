import math

import numpy as np
import nurbforge as nf
import pytest


def test_version():
    assert nf.__version__ == "0.1.0"


def test_circle_exact():
    c = nf.make_circle((1, 2, 3), (0, 0, 1), 2.0)
    pts = c.sample(200)
    r = np.linalg.norm(pts - np.array([1, 2, 3]), axis=1)
    assert np.allclose(r, 2.0, atol=1e-12)


def test_curve_from_numpy_and_knot_insertion():
    ctrl = np.array([[0, 0, 0], [1, 2, 0], [3, 2, 1], [4, 0, 0]], dtype=float)
    c = nf.NurbsCurve.from_points(3, ctrl)
    c2 = c.insert_knot(0.4, 2)
    assert len(c2.control_points) == 6
    for u in np.linspace(0, 1, 11):
        assert nf.distance(c(u), c2(u)) < 1e-12


def test_invalid_input_raises():
    with pytest.raises(nf.GeometryError):
        nf.NurbsCurve(2, [0, 0, 1, 1], np.zeros((3, 3)))
    with pytest.raises(ValueError):
        nf.make_circle((0, 0, 0), (0, 0, 0), 1.0)


def test_sphere_by_revolution():
    arc = nf.make_arc((0, 0, 0), (0, 0, -1), (1, 0, 0), 2.0, 0.0, math.pi)
    s = nf.revolve(arc, (0, 0, 0), (0, 0, 1), 2 * math.pi)
    grid = s.sample(15, 15)
    assert grid.shape == (15, 15, 3)
    assert np.allclose(np.linalg.norm(grid, axis=2), 2.0)


def test_mirror_and_polar_array():
    line = nf.make_line((1, 0, 0), (2, 0, 0))
    m = line.transformed(nf.Transform.mirror((0, 0, 0), (1, 0, 0)))
    assert m(0.0).x == pytest.approx(-1.0)
    copies = nf.apply_pattern(line, nf.polar_pattern((0, 0, 0), (0, 0, 1), 6, 2 * math.pi))
    assert len(copies) == 6
    assert copies[3](0.0).x == pytest.approx(-1.0)


def test_surface_normal_and_iso_curve():
    s = nf.make_rectangle((0, 0, 0), (1, 0, 0), (0, 1, 0))
    n = s.normal(0.5, 0.5)
    assert tuple(n) == pytest.approx((0, 0, 1))
    assert nf.distance(s.iso_curve_u(0.3)(0.6), s(0.3, 0.6)) < 1e-12


def test_curve_refinement_and_bezier_decomposition():
    c = nf.make_circle((0, 0, 0), (0, 0, 1), 2.0)
    r = c.refine_knots([0.6, 0.1, 0.6])
    assert len(r.control_points) == len(c.control_points) + 3
    for u in np.linspace(0, 1, 51):
        assert nf.distance(c(u), r(u)) < 1e-12

    segs = c.decompose_bezier()
    assert len(segs) == 4
    for s in segs:
        a, b = s.domain
        for u in np.linspace(a, b, 7):
            assert nf.distance(s(u), c(u)) < 1e-12

    with pytest.raises(nf.GeometryError):
        c.refine_knots([1.5])


def test_surface_refinement_and_bezier_decomposition():
    meridian = nf.make_arc((0, 0, 0), (0, 0, -1), (1, 0, 0), 1.0, 0.0, math.pi)
    sphere = nf.revolve(meridian, (0, 0, 0), (0, 0, 1), 2 * math.pi)
    r = sphere.refine_knots_u([0.1, 0.9]).refine_knots_v([0.3])
    assert (r.num_u, r.num_v) == (sphere.num_u + 2, sphere.num_v + 1)
    assert nf.distance(r(0.37, 0.61), sphere(0.37, 0.61)) < 1e-12

    patches = sphere.decompose_bezier()
    assert len(patches) == 4 * 2
    for p in patches:
        (u0, u1), (v0, v1) = p.domain_u, p.domain_v
        u, v = 0.5 * (u0 + u1), 0.5 * (v0 + v1)
        assert nf.distance(p(u, v), sphere(u, v)) < 1e-12


def test_degree_elevation():
    c = nf.make_circle((0, 0, 0), (0, 0, 1), 1.5)
    e = c.elevate_degree(2)
    assert e.degree == 4 and e.is_clamped
    for u in np.linspace(0, 1, 41):
        assert nf.distance(c(u), e(u)) < 1e-12

    meridian = nf.make_arc((0, 0, 0), (0, 0, -1), (1, 0, 0), 1.0, 0.0, math.pi)
    sphere = nf.revolve(meridian, (0, 0, 0), (0, 0, 1), 2 * math.pi)
    s = sphere.elevate_degree_u().elevate_degree_v(2)
    assert (s.degree_u, s.degree_v) == (3, 4)
    assert nf.distance(s(0.3, 0.7), sphere(0.3, 0.7)) < 1e-12


def test_ruled_surface_between_different_degrees():
    line = nf.make_line((-1, 0, 0), (1, 0, 0))
    arc = nf.make_arc((0, 0, 2), (-1, 0, 0), (0, 1, 0), 1.0, 0.0, math.pi)
    s = nf.ruled(line, arc)
    for u in np.linspace(0, 1, 11):
        assert nf.distance(s(u, 0.0), line(u)) < 1e-12
        assert nf.distance(s(u, 1.0), arc(u)) < 1e-12


def test_knot_removal_and_degree_reduction():
    c = nf.make_circle((0, 0, 0), (0, 0, 1), 2.0)
    r = c.insert_knot(0.6, 2).remove_knot(0.6, 2)
    assert np.allclose(r.control_points, c.control_points, atol=1e-12)

    same, removed = c.try_remove_knot(0.25)
    assert len(same.control_points) == len(c.control_points)
    assert removed == 0
    with pytest.raises(nf.GeometryError):
        c.remove_knot(0.25)

    back = c.elevate_degree(1).reduce_degree()
    assert back.degree == 2
    assert np.allclose(back.control_points, c.control_points, atol=1e-12)
    assert np.allclose(back.knots, c.knots, atol=1e-12)

    wavy = nf.NurbsCurve.from_points(3, np.array(
        [[0, 0, 0], [1, 3, 0], [2, -2, 0], [3, 3, 1], [4, -1, 0]], dtype=float))
    with pytest.raises(nf.GeometryError):
        wavy.reduce_degree()

    meridian = nf.make_arc((0, 0, 0), (0, 0, -1), (1, 0, 0), 1.0, 0.0, math.pi)
    sphere = nf.revolve(meridian, (0, 0, 0), (0, 0, 1), 2 * math.pi)
    s = sphere.elevate_degree_v().reduce_degree_v().refine_knots_u([0.1]).remove_knot_u(0.1)
    assert (s.degree_v, s.num_u, s.num_v) == (2, sphere.num_u, sphere.num_v)
    assert nf.distance(s(0.3, 0.7), sphere(0.3, 0.7)) < 1e-12


def test_split_and_join():
    c = nf.make_circle((0, 0, 0), (0, 0, 1), 2.0)
    left, right = c.split(0.4)
    assert left.domain == (0.0, 0.4) and right.domain == (0.4, 1.0)
    back = nf.join(left, right)
    assert np.allclose(back.control_points, c.control_points, atol=1e-12)
    assert np.allclose(back.knots, c.knots, atol=1e-12)

    s = c.sub_curve(0.1, 0.3)
    assert nf.distance(s(0.2), c(0.2)) < 1e-12

    quarter = c.sub_curve(0.0, 0.25)
    line = nf.make_line(tuple(c(0.25)), (5, 5, 0))
    j = nf.join(quarter, line)
    assert j.degree == 2 and j.domain == (0.0, 1.25)
    assert nf.distance(j(0.1), c(0.1)) < 1e-12
    assert nf.distance(j(1.25), nf.Vec3(5, 5, 0)) < 1e-12
    with pytest.raises(nf.GeometryError):
        nf.join(line, quarter)  # line ends at (5,5,0), the arc does not start there

    meridian = nf.make_arc((0, 0, 0), (0, 0, -1), (1, 0, 0), 1.0, 0.0, math.pi)
    sphere = nf.revolve(meridian, (0, 0, 0), (0, 0, 1), 2 * math.pi)
    a, b = sphere.split_u(0.3)
    ju = nf.join_u(a, b)
    assert ju.num_u == sphere.num_u
    assert nf.distance(ju(0.6, 0.4), sphere(0.6, 0.4)) < 1e-12
    patch = sphere.sub_surface(0.1, 0.2, 0.3, 0.6)
    assert nf.distance(patch(0.15, 0.5), sphere(0.15, 0.5)) < 1e-12


def test_surface_derivatives_and_pole_normals():
    meridian = nf.make_arc((0, 0, 0), (0, 0, -1), (1, 0, 0), 2.0, 0.0, math.pi)
    sphere = nf.revolve(meridian, (0, 0, 0), (0, 0, 1), 2 * math.pi)
    D = sphere.derivatives(0.3, 0.4, 2)
    assert len(D) == 3 and len(D[0]) == 3
    h = 1e-6
    fd = (np.array(sphere(0.3 + h, 0.4)) - np.array(sphere(0.3 - h, 0.4))) / (2 * h)
    assert np.allclose(fd, np.array(D[1][0]), atol=1e-6)

    sign = np.sign(np.dot(np.array(sphere.normal(0.3, 0.5)), np.array(sphere(0.3, 0.5))))
    for v in (0.0, 1.0):  # poles
        n = np.array(sphere.normal(0.25, v))
        assert np.allclose(n, sign * np.array(sphere(0.25, v)) / 2.0, atol=1e-9)


def test_closest_point():
    c = nf.make_circle((0, 0, 0), (0, 0, 1), 2.0)
    r = nf.closest_point(c, (5, 0, 1))
    assert isinstance(r, nf.CurvePoint)
    assert abs(r.distance - math.hypot(3, 1)) < 1e-12
    assert nf.distance(r.point, nf.Vec3(2, 0, 0)) < 1e-9
    assert nf.distance(c(r.u), r.point) < 1e-14

    meridian = nf.make_arc((0, 0, 0), (0, 0, -1), (1, 0, 0), 2.0, 0.0, math.pi)
    sphere = nf.revolve(meridian, (0, 0, 0), (0, 0, 1), 2 * math.pi)
    s = nf.closest_point(sphere, (0, 0, 7), nf.Tolerance(length=1e-12))
    assert abs(s.distance - 5.0) < 1e-12  # pole
    q = sphere(0.3, 0.6)
    loc = nf.closest_point_local(sphere, q, 0.28, 0.62)
    assert loc.distance < 1e-10


def test_properties():
    c = nf.make_circle((0, 0, 0), (0, 0, 1), 2.0)
    assert abs(nf.length(c) - 4 * math.pi) < 1e-11
    assert abs(nf.curvature(c, 0.3) - 0.5) < 1e-12
    assert abs(nf.arc_length(c, 0.0, 0.25) - math.pi) < 1e-11
    u = nf.parameter_at_length(c, 1.0)
    assert abs(nf.arc_length(c, 0.0, u) - 1.0) < 1e-9

    meridian = nf.make_arc((0, 0, 0), (0, 0, -1), (1, 0, 0), 2.0, 0.0, math.pi)
    sphere = nf.revolve(meridian, (0, 0, 0), (0, 0, 1), 2 * math.pi)
    assert abs(nf.area(sphere) - 16 * math.pi) < 1e-9
    k = nf.surface_curvature(sphere, 0.2, 0.4)
    assert abs(k.gaussian - 0.25) < 1e-12 and abs(abs(k.mean) - 0.5) < 1e-12
