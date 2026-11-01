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


def test_fitting():
    t = np.linspace(0, 1, 15)
    pts = np.stack([t, t**3, 1 - t**2], axis=1)
    c = nf.interpolate_curve(pts, 3, nf.Parametrization.Uniform)
    for tt in (0.13, 0.5, 0.91):
        assert nf.distance(c(tt), nf.Vec3(tt, tt**3, 1 - tt**2)) < 1e-12

    d = nf.interpolate_curve(pts, 3, end_derivatives=((1, 0, 0), (1, 3, -2)))
    assert nf.distance(d.derivatives(0.0, 1)[1], nf.Vec3(1, 0, 0)) < 1e-10

    a = nf.approximate_curve(pts, 3, 6)
    assert len(a.control_points) == 6
    f = nf.approximate_curve_tol(pts, 3, 1e-3)
    ub = nf.data_parameters(pts)
    assert max(nf.distance(f(u), nf.Vec3(*p)) for u, p in zip(ub, pts)) <= 1e-3

    x, y = np.meshgrid(np.linspace(0, 1, 8), np.linspace(0, 1, 6), indexing="ij")
    grid = np.stack([x, y, np.sin(x) * np.cos(y)], axis=2)
    s = nf.interpolate_surface(grid, 3, 3, nf.Parametrization.Uniform)
    assert nf.distance(s(x[3, 2], y[3, 2]), nf.Vec3(*grid[3, 2])) < 1e-12
    w = nf.approximate_surface(grid, 2, 2, 4, 4)
    assert (w.num_u, w.num_v) == (4, 4)


def test_extend():
    arc = nf.make_arc((0, 0, 0), (1, 0, 0), (0, 1, 0), 1.0, 0.0, math.pi / 2)
    e = nf.extend(arc, nf.CurveEnd.End, math.pi / 4)
    assert abs(nf.length(e) - 3 * math.pi / 4) < 1e-9
    assert abs(nf.norm(e(e.domain[1])) - 1.0) < 1e-12  # still on the circle
    line = nf.extend(nf.make_line((0, 0, 0), (1, 0, 0)), nf.CurveEnd.Start, 2.0, nf.Extension.Linear)
    assert nf.distance(line(line.domain[0]), nf.Vec3(-2, 0, 0)) < 1e-12

    cyl = nf.extrude(nf.make_circle((0, 0, 0), (0, 0, 1), 1.0), (0, 0, 2))
    s = nf.extend(cyl, nf.SurfaceSide.VMax, 1.0)
    assert abs(s(0.2, s.domain_v[1]).z - 3.0) < 1e-12


def test_loft():
    secs = [nf.make_circle((0, 0, z), (0, 0, 1), 1.0 + 0.5 * z) for z in (0.0, 1.0, 2.0, 3.0)]
    s = nf.loft(secs)
    p = s(0.37, 0.61)
    assert abs(math.hypot(p.x, p.y) - (1.0 + 0.5 * p.z)) < 1e-12  # exact cone
    v = nf.loft_parameters(secs)
    assert v[0] == 0.0 and v[-1] == 1.0
    c = nf.make_compatible([nf.make_line((0, 0, 0), (1, 0, 0)), secs[0]])
    assert c[0].degree == c[1].degree and list(c[0].knots) == list(c[1].knots)


def test_sweep():
    profile = nf.make_circle((5, 0, 0), (0, 1, 0), 1.0)
    path = nf.make_arc((0, 0, 0), (1, 0, 0), (0, 1, 0), 5.0, 0.0, math.pi)
    s = nf.sweep(profile, path, nf.SweepOptions(tolerance=1e-7))
    p = s(0.3, 0.6)
    assert abs(math.hypot(math.hypot(p.x, p.y) - 5.0, p.z) - 1.0) < 1e-6  # torus
    t = nf.sweep_translational(nf.make_line((0, 0, 0), (1, 0, 0)), nf.make_line((0, 0, 0), (0, 2, 0)))
    assert nf.distance(t(1.0, 1.0), nf.Vec3(1, 2, 0)) < 1e-14
    r = nf.sweep2(nf.make_line((0, 0, 0), (1, 0, 0)), nf.make_line((0, 0, 0), (0, 5, 0)),
                  nf.make_line((1, 0, 0), (2, 5, 0)))
    assert nf.distance(r(1.0, 1.0), nf.Vec3(2, 5, 0)) < 1e-6


def test_coons_gordon():
    p00, p10, p01, p11 = (0, 0, 0), (2, 0, 1), (0, 3, 0), (2, 3, 4)
    s = nf.coons(nf.make_line(p00, p10), nf.make_line(p01, p11),
                 nf.make_line(p00, p01), nf.make_line(p10, p11))
    b = nf.make_bilinear(p00, p10, p01, p11)
    assert nf.distance(s(0.3, 0.8), b(0.3, 0.8)) < 1e-12
    uc = [b.iso_curve_v(v) for v in (0.0, 0.5, 1.0)]
    vc = [b.iso_curve_u(u) for u in (0.0, 0.4, 1.0)]
    g = nf.gordon(uc, vc)
    assert nf.distance(g(0.7, 0.2), b(0.7, 0.2)) < 1e-10


def test_optics():
    sphere = nf.Asphere(curvature=1 / 50.0)
    s = nf.asphere_surface(sphere, 20.0)
    p = s(0.3, 0.7)
    assert abs(math.sqrt(p.x**2 + p.y**2 + (p.z - 50.0) ** 2) - 50.0) < 1e-12  # exact
    asph = nf.Asphere(curvature=1 / 40.0, conic=-0.8, even=[1e-6, -2e-9])
    a = nf.asphere_surface(asph, 15.0, tolerance=1e-8)
    q = a(0.1, 0.55)
    assert abs(q.z - asph.sag(math.hypot(q.x, q.y))) < 2e-8

    assert abs(nf.zernike(4, nf.ZernikeIndexing.Fringe, 0.5, 0.0) - (2 * 0.25 - 1)) < 1e-15
    assert nf.zernike_nm(11, nf.ZernikeIndexing.Noll) == (4, 0)
    z = nf.ZernikeSag(base=sphere, norm_radius=10.0, coeffs=[0, 0.01, 0, 0.02])
    zs = nf.zernike_surface(z, tolerance=1e-8)
    assert abs(zs(3.0, -4.0).z - z.sag(3.0, -4.0)) < 2e-8
    xy = nf.XYPolynomialSag(norm_radius=5.0, terms=[(2, 0, 0.1), (1, 2, 0.05)])
    assert abs(nf.xy_polynomial_surface(xy)(1.0, 2.0).z - xy.sag(1.0, 2.0)) < 1e-13
    g = nf.sag_surface(lambda x, y: math.sin(x) * math.cos(y), 0, 1, 0, 1, tolerance=1e-9)
    assert abs(g(0.3, 0.6).z - math.sin(0.3) * math.cos(0.6)) < 2e-9


def test_intersections():
    circle = nf.make_circle((0, 0, 0), (0, 0, 1), 2.0)
    hits = nf.intersect(circle, nf.make_line((-3, 1, 0), (3, 1, 0)))
    assert len(hits) == 2 and all(abs(h.point.y - 1) < 1e-12 for h in hits)
    meridian = nf.make_arc((0, 0, 0), (0, 0, -1), (1, 0, 0), 2.0, 0.0, math.pi)
    sphere = nf.revolve(meridian, (0, 0, 0), (0, 0, 1), 2 * math.pi)
    sh = nf.intersect(nf.make_line((0, 0, -3), (0, 0, 3)), sphere)
    assert [round(h.point.z, 12) for h in sh] == [-2.0, 2.0]


def test_surface_surface_intersection():
    meridian = nf.make_arc((0, 0, 0), (0, 0, -1), (1, 0, 0), 2.0, 0.0, math.pi)
    sphere = nf.revolve(meridian, (0, 0, 0), (0, 0, 1), 2 * math.pi)
    plane = nf.make_rectangle((-3, -3, 0.8), (6, 0, 0), (0, 6, 0))
    res = nf.intersect(sphere, plane)
    total = sum(nf.length(c.curve) for c in res.curves)
    assert abs(total - 2 * math.pi * math.sqrt(4 - 0.64)) < 1e-5
    c = res.curves[0]
    w = c.pcurve1(0.5 * sum(c.curve.domain))
    assert abs(sphere(w.x, w.y).z - 0.8) < 1e-6
    assert res.points == []
