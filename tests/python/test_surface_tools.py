import math

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
