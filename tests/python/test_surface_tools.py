import math

import numpy as np
import nurbforge as nf
import pytest


def test_rebuild():
    circle = nf.make_circle((0, 0, 0), (0, 0, 1), 1.0)
    a, b = nf.rebuild_curve(circle, 12), nf.rebuild_curve(circle, 40)
    assert b.max_error < a.max_error < 1e-2
    assert len(b.curve.control_points) == 40
    patch = nf.make_bilinear((0, 0, 0), (2, 0, 1), (0, 1, 0), (2, 1, 3))
    r = nf.rebuild_surface(patch, 5, 6)
    assert r.max_error < 1e-12
    with pytest.raises(nf.GeometryError):
        nf.rebuild_curve(circle, 2)


def test_tapered_extrusion_ribbon_cone():
    square = nf.make_rectangle_curve((0, 0, 0), (2, 0, 0), (0, 2, 0))
    ang = math.radians(10)
    faces = nf.extrude_tapered(square, (0, 0, 1), ang)
    assert len(faces) == 4
    solid = nf.cap(nf.Shell.sew([nf.TrimmedSurface(f) for f in faces]))
    a1, a2 = 4.0, (2 + 2 * math.tan(ang)) ** 2
    assert abs(solid.volume()) == pytest.approx((a1 + a2 + math.sqrt(a1 * a2)) / 3, abs=1e-9)
    strip = nf.ribbon(nf.make_polyline([(0, 0, 0), (2, 0, 0), (2, 1, 0)]), 0.5)
    assert sum(nf.TrimmedSurface(s).area() for s in strip) == pytest.approx(1.25, abs=1e-9)
    cone = nf.extrude_to_point(nf.make_circle((0, 0, 0), (0, 0, 1), 1.0), (0, 0, 2))
    assert nf.TrimmedSurface(cone).area() == pytest.approx(math.pi * math.sqrt(5), abs=1e-8)


def test_coons_exact_and_triangle():
    a0 = nf.make_arc((0, 0, 0), (1, 0, 0), (0, 1, 0), 1.0, 0, math.pi / 2)
    a1 = nf.make_arc((0, 0, 1), (1, 0, 0), (0, 1, 0), 1.0, 0, math.pi / 2)
    s = nf.coons(a0, a1, nf.make_line((1, 0, 0), (1, 0, 1)), nf.make_line((0, 1, 0), (0, 1, 1)))
    for i in range(6):
        for j in range(6):
            p = s(i / 5, j / 5)
            assert math.hypot(p.x, p.y) == pytest.approx(1, abs=1e-12)
    arcs = [nf.make_arc((0, 0, 0), x, y, 1.0, 0, math.pi / 2)
            for x, y in (((1, 0, 0), (0, 1, 0)), ((0, 1, 0), (0, 0, 1)), ((0, 0, 1), (1, 0, 0)))]
    t = nf.triangle_patch(*arcs)
    p = t(0.5, 1.0)
    assert (p.x, p.y, p.z) == pytest.approx((0, 0, 1), abs=1e-12)


def test_multi_profile_sweep_and_rail_revolve():
    path = nf.make_line((0, 0, 0), (0, 0, 4))
    prof = [nf.make_circle((0, 0, 0), (0, 0, 1), 1.0), nf.make_circle((0, 0, 4), (0, 0, 1), 2.0)]
    s = nf.sweep(prof, path, options=nf.SweepOptions(tolerance=1e-8))
    p = s(0.3, 0.5)
    assert math.hypot(p.x, p.y) == pytest.approx(1.5, abs=1e-7)
    rail = nf.make_circle((0, 0, 0), (0, 0, 1), 1.0).transformed(
        nf.Transform.scaling((0, 0, 0), (2, 1, 1)))
    disc = nf.rail_revolve(nf.make_line((0, 0, 0), (2, 0, 0)), rail, (0, 0, 0), (0, 0, 1),
                           options=nf.SweepOptions(tolerance=1e-8))
    assert nf.TrimmedSurface(disc).area() == pytest.approx(2 * math.pi, abs=1e-6)


def _plane_strip(x0, z0, slope):
    pts = [[(x0 + i / 3, j / 3, z0 + slope * (x0 + i / 3) ** 2 + 0.1 * j / 3) for j in range(4)]
           for i in range(4)]
    return nf.interpolate_surface(pts)


def test_blend_and_match_surfaces():
    a, b = _plane_strip(0, 0, 0.2), _plane_strip(2, 1, -0.1)
    for c, deg in ((0, 1), (1, 3), (2, 5)):
        bl = nf.blend_surfaces(a, nf.SurfaceSide.UMax, b, nf.SurfaceSide.UMin, continuity=c)
        assert bl.degree_v == deg
        for t in (0.0, 0.3, 1.0):
            assert nf.distance(bl(t, 0.0), a(1.0, t)) < 1e-12
            assert nf.distance(bl(t, 1.0), b(0.0, t)) < 1e-12
            if c >= 1:
                assert nf.norm(nf.cross(bl.normal(t, 0.0), a.normal(1.0, t))) < 1e-10
                assert nf.norm(nf.cross(bl.normal(t, 1.0), b.normal(0.0, t))) < 1e-10
            if c == 2:
                k0 = nf.surface_curvature(bl, t, 0.0).mean
                k1 = nf.surface_curvature(a, 1.0, t).mean
                assert abs(abs(k0) - abs(k1)) < 1e-8
    m = nf.match_surface(b, nf.SurfaceSide.UMin, a, nf.SurfaceSide.UMax, continuity=1)
    for t in (0.0, 0.5, 1.0):
        assert nf.distance(m(0.0, t), a(1.0, t)) < 1e-12
        assert nf.norm(nf.cross(m.normal(0.0, t), a.normal(1.0, t))) < 1e-10
        assert nf.distance(m(1.0, t), b(1.0, t)) < 1e-12
    with pytest.raises(nf.GeometryError):
        nf.blend_surfaces(a, nf.SurfaceSide.UMax, b, nf.SurfaceSide.UMin, continuity=3)


def test_fit_surface_scattered():
    rng = np.random.default_rng(3)
    xy = rng.random((600, 2)) * [3, 2]
    pts = np.column_stack([xy, 0.1 * np.sin(xy[:, 0]) * np.cos(xy[:, 1])])
    f = nf.fit_surface(pts)
    assert f.max_error < 1e-4
    assert f.parameters.shape == (600, 2)
    assert "FittedSurface" in repr(f)
    edges = [nf.make_line((0, 0, 0), (1, 0, 0)), nf.make_line((0, 1, 0), (1, 1, 0)),
             nf.make_line((0, 0, 0), (0, 1, 0)), nf.make_line((1, 0, 0), (1, 1, 0))]
    uv = rng.random((400, 2))
    bump = np.column_stack([uv, 0.1 * np.sin(np.pi * uv[:, 0]) * np.sin(np.pi * uv[:, 1])])
    g = nf.fit_surface(bump, boundary=edges, num_u=8, num_v=8)
    assert abs(g.surface(0.5, 0.0).z) < 1e-15
    with pytest.raises(nf.GeometryError):
        nf.fit_surface(pts[:2])


def _xyz(p):
    return (p.x, p.y, p.z)


def test_sweep2_profiles_and_scaling():
    rail1, rail2 = nf.make_line((0, 0, 0), (0, 10, 0)), nf.make_line((2, 0, 0), (4, 10, 0))
    arch = nf.make_arc((1, 0, 0), (-1, 0, 0), (0, 0, 1), 1.0, 0.0, math.pi)
    uni = nf.sweep2(arch, rail1, rail2)
    non = nf.sweep2([arch], rail1, rail2, scaling=nf.Sweep2Scaling.NonUniform)
    assert np.allclose(_xyz(uni(0.5, 1.0)), (2, 10, 2), atol=1e-9)  # scaled by 2 at the rail ends
    assert np.allclose(_xyz(non(0.5, 1.0)), (2, 10, 1), atol=1e-9)  # height kept
    # a second, flat profile at the rail ends: the arch fades out
    flat = nf.make_line((0, 10, 0), (4, 10, 0))
    s = nf.sweep2([arch, flat], rail1, rail2)
    assert np.allclose(_xyz(s(0.5, 0.0)), (1, 0, 1), atol=1e-9)
    assert abs(s(0.5, 1.0).z) < 1e-9
    assert 0.0 < s(0.5, 0.5).z < 1.5
    with pytest.raises(nf.GeometryError):
        nf.sweep2([arch, flat], rail1, rail2, params=[0.3, 0.3])


def test_variable_offset():
    plane = nf.make_bilinear((0, 0, 0), (2, 0, 0), (0, 3, 0), (2, 3, 0))
    wedge = nf.offset_surface_variable(plane, [0.1, 0.5, 0.1, 0.5])
    du, dv = plane.domain_u, plane.domain_v
    um = 0.5 * (du[0] + du[1])
    assert abs(abs(wedge(um, dv[0]).z) - 0.3) < 1e-9
    bump = nf.offset_surface_variable(
        plane, lambda u, v: 0.2 * math.sin(math.pi * (u - du[0]) / (du[1] - du[0])), 1e-8)
    assert abs(abs(bump(um, dv[1]).z) - 0.2) < 1e-8
    assert abs(bump(du[0], dv[0]).z) < 1e-8
    # on a unit sphere a distance of 5 is possible outwards only
    sphere = nf.revolve(nf.make_arc((0, 0, 0), (1, 0, 0), (0, 0, 1), 1.0, -math.pi / 2,
                                    math.pi / 2), (0, 0, 0), (0, 0, 1), 2 * math.pi)
    refused = 0
    for d in (5.0, -5.0):
        try:
            nf.offset_surface_variable(sphere, lambda u, v, d=d: d, 1e-5)
        except nf.GeometryError:
            refused += 1
    assert refused == 1
