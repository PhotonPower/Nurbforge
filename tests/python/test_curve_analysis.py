"""Curve analysis: curvature, curvature graph, direction, continuity, duplicates."""
import math

import nurbforge as nf
import pytest


def test_length_degree_curvature():
    circle = nf.make_circle((1, 2, 3), (0, 0, 1), 2.5)
    assert nf.length(circle) == pytest.approx(5 * math.pi, abs=1e-9)
    assert circle.degree == 2
    k = nf.curve_curvature(circle, 0.3)
    assert k.curvature == pytest.approx(0.4, abs=1e-12)
    assert k.radius == pytest.approx(2.5, abs=1e-11)
    assert nf.distance(k.center, nf.Vec3(1, 2, 3)) < 1e-11
    assert abs(nf.dot(k.tangent, k.normal)) < 1e-12
    flat = nf.curve_curvature(nf.make_line((0, 0, 0), (1, 1, 0)), 0.5)
    assert flat.curvature == 0.0 and math.isinf(flat.radius)


def test_curvature_graph():
    g = nf.curvature_graph(nf.make_circle((0, 0, 0), (0, 0, 1), 2.0), scale=3.0,
                           samples_per_span=8)
    assert len(g.points) == len(g.tips) == len(g.parameters) == len(g.curvatures)
    assert g.max_curvature == pytest.approx(0.5, abs=1e-12)
    assert all(abs(nf.norm(t) - 3.5) < 1e-11 for t in g.tips)


def test_direction():
    ccw = nf.make_circle((0, 0, 0), (0, 0, 1), 2.0)
    d = nf.curve_direction(ccw)
    assert d.closed and d.planar and d.orientation == 1
    assert d.signed_area == pytest.approx(4 * math.pi, abs=1e-9)
    assert nf.curve_direction(ccw.reversed()).orientation == -1
    assert nf.curve_direction(ccw, plane_normal=(0, 0, -1)).orientation == -1
    arc = nf.make_arc((0, 0, 0), (1, 0, 0), (0, 1, 0), 1.0, 0.0, math.pi / 2)
    o = nf.curve_direction(arc)
    assert not o.closed and o.orientation == 0
    assert nf.distance(o.start_tangent, nf.Vec3(0, 1, 0)) < 1e-12
    assert nf.is_closed(ccw) and not nf.is_closed(arc) and nf.is_planar(arc)


def test_continuity():
    line = nf.make_line((-2, 1, 0), (0, 1, 0))
    arc = nf.make_arc((0, 0, 0), (0, 1, 0), (1, 0, 0), 1.0, 0.0, math.pi / 2)
    k = nf.continuity(line, arc)
    assert k.geometric == 1 and k.gap < 1e-12 and k.tangent_angle < 1e-12
    assert k.curvature_a == 0.0 and k.curvature_b == pytest.approx(1.0)
    assert nf.continuity(line, nf.CurveEnd.End, arc, nf.CurveEnd.Start).geometric == 1
    corner = nf.continuity(line, nf.make_line((0, 1, 0), (0, 3, 0)))
    assert corner.geometric == 0 and corner.tangent_angle == pytest.approx(math.pi / 2)
    assert nf.continuity(line, nf.make_line((0, 2, 0), (0, 3, 0))).geometric == -1
    breaks = nf.continuity_breaks(nf.join(line, arc))
    assert len(breaks) == 1 and breaks[0].continuity.geometric == 1
    assert not nf.continuity_breaks(nf.make_circle((0, 0, 0), (0, 0, 1), 1.0))
    assert len(nf.continuity_breaks(nf.make_rectangle_curve((0, 0, 0), (2, 0, 0), (0, 1, 0)))) == 4


def test_duplicates():
    line = nf.make_line((0, 0, 0), (2, 1, 0))
    circle = nf.make_circle((0, 0, 0), (0, 0, 1), 1.0)
    curves = [circle, line, line.reversed(), nf.make_line((0, 0, 0), (1, 0.5, 0)),
              circle.elevate_degree()]
    assert nf.find_duplicate_curves(curves, 1e-9) == [[0, 4], [1, 2]]
    assert nf.find_duplicate_curves([circle, line]) == []
