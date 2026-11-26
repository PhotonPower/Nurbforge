"""Silhouettes, hidden-line removal, Make2D and outlines."""
import math

import nurbforge as nf
import pytest


def _length(curves):
    return sum(nf.length(c.curve) for c in curves)


def test_silhouette_of_a_sphere():
    sphere = nf.make_sphere((0, 0, 0), 2.0)
    d = nf.Vec3(1, 2, 3)
    curves = nf.silhouette(sphere, direction=d)
    assert _length(curves) == pytest.approx(4 * math.pi, abs=1e-5)
    for c in curves:
        assert c.face == 0 and not c.crease and not c.cusps
        p = c.curve(0.5 * (c.curve.domain[0] + c.curve.domain[1]))
        assert abs(nf.norm(p) - 2.0) < 1e-6
        assert abs(nf.dot(p, d)) < 1e-5
        q = c.pcurve(c.pcurve.domain[0])
        assert nf.distance(sphere.faces[0].surface(q.x, q.y), c.curve(c.curve.domain[0])) < 1e-6
    # from an eye point at distance 5: a small circle
    eye_curves = nf.silhouette(sphere.faces[0].surface, eye=(5, 0, 0))
    assert _length(eye_curves) == pytest.approx(4 * math.pi * math.sqrt(1 - 4 / 25), abs=1e-5)
    with pytest.raises(nf.GeometryError):
        nf.silhouette(sphere)
    with pytest.raises(nf.GeometryError):
        nf.silhouette(sphere, direction=(1, 0, 0), eye=(5, 0, 0))


def test_visibility_behind_a_box():
    box = nf.make_box((0, 0, 0), (1, 1, 1))
    line = nf.make_line((-1, 0.5, -1), (2, 0.5, -1))
    pieces = nf.visibility([line], [box], direction=(0, 0, -1))
    assert [p.visible for p in pieces] == [True, False, True]
    assert pieces[1].t0 == pytest.approx(1 / 3, abs=1e-7)
    assert pieces[1].t1 == pytest.approx(2 / 3, abs=1e-7)
    assert all(p.source == 0 for p in pieces)


def test_make2d_and_outline_of_a_plate_with_a_hole():
    plate = nf.boolean(nf.make_box((0, 0, 0), (4, 2, 1)),
                       nf.make_cylinder((2, 1, -1), (0, 0, 1), 0.5, 3), nf.BooleanOp.Difference)
    top = nf.make2d([plate], (0, 0, -1), s_axis=(1, 0, 0))
    assert not top.hidden
    assert _length(top.visible) == pytest.approx(12 + math.pi, abs=1e-5)
    assert nf.distance(top.t_axis, nf.Vec3(0, 1, 0)) < 1e-12
    front = nf.make2d([plate], (0, 1, 0), s_axis=(1, 0, 0))
    assert _length(front.visible) == pytest.approx(10.0, abs=1e-5)
    assert _length(front.hidden) == pytest.approx(2.0, abs=1e-5)  # the hole
    assert all(c.kind == nf.DrawingCurveKind.Silhouette for c in front.hidden)
    assert not nf.make2d([plate], (0, 1, 0), hidden=False).hidden
    contour = nf.outline([plate], (0, 0, -1))
    assert len(contour.loops) == 2
    assert nf.region_area(contour.loops) == pytest.approx(8 - math.pi / 4, abs=1e-6)


def test_isometric_cube():
    cube = nf.make_box((0, 0, 0), (1, 1, 1))
    d = nf.make2d([cube], (-1, -1, -1))
    assert len(d.visible) == 9 and len(d.hidden) == 3
    hexagon = nf.outline([cube], (-1, -1, -1))
    assert len(hexagon.loops) == 1 and len(hexagon.loops[0]) == 6
    assert nf.region_area(hexagon.loops) == pytest.approx(math.sqrt(3), abs=1e-9)
