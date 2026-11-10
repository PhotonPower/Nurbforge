"""Hatch and contour lines of section loops from Python."""

import math

import nurbforge as nf
import pytest


def box(h: float) -> nf.Shell:
    s = 2 * h
    rect = nf.make_rectangle
    return nf.Shell.sew([
        nf.TrimmedSurface(rect((-h, -h, -h), (s, 0, 0), (0, s, 0))),
        nf.TrimmedSurface(rect((-h, -h, -h), (s, 0, 0), (0, 0, s))),
        nf.TrimmedSurface(rect((-h, h, -h), (s, 0, 0), (0, 0, s))),
        nf.TrimmedSurface(rect((-h, -h, -h), (0, s, 0), (0, 0, s))),
        nf.TrimmedSurface(rect((h, -h, -h), (0, s, 0), (0, 0, s))),
        nf.TrimmedSurface(rect((-h, -h, h), (s, 0, 0), (0, s, 0))),
    ])


def test_circle_chords_and_area():
    r = 2.0
    circle = nf.make_circle((0, 0, 0), (0, 0, 1), r)
    h = 0.05
    segments, contours = nf.hatch([circle], (1, 0), h, phase=0.01)
    assert contours == []
    for s in segments:
        y = s.start[1]
        assert math.isclose(y, 0.01 + s.line * h, abs_tol=1e-14)
        assert math.isclose(s.length, 2 * math.sqrt(r * r - y * y), abs_tol=1e-12)
    total = sum(s.length for s in segments)
    assert abs(h * total - math.pi * r * r) < 2 * math.pi * r * h


def test_section_of_a_box_with_contours():
    sec = nf.section(box(1.0), nf.Plane((0, 0, 0.25), (0, 0, 1)))
    # the PlanarSection itself or its loops
    for loops in (sec, sec.loops):
        segments, contours = nf.hatch(loops, (1, 1), 0.1, contours=[0.1, 0.2])
        assert all(s.length > 0 and s.region == 0 for s in segments)
        assert [(c.loop, c.index, c.distance) for c in contours] == [(0, 0, 0.1), (0, 1, 0.2)]
        assert len(contours[1].curves) == 4


def test_ring_regions_and_tangent_lines():
    outer = nf.make_circle((0, 0, 0), (0, 0, 1), 2.0)
    hole = nf.make_circle((0, 0, 0), (0, 0, 1), 1.0).reversed()
    segments, _ = nf.hatch([[outer], [hole]], (1, 0, 0), 0.5)
    # y = +-1 touch the hole, y = +-2 touch the outer circle: no zero-length segments
    assert all(s.length > 0.1 for s in segments)
    assert {s.line for s in segments} == {-3, -2, -1, 0, 1, 2, 3}
    assert sum(1 for s in segments if s.line == 0) == 2
    assert {s.region for s in segments} == {0}


def test_errors():
    circle = nf.make_circle((0, 0, 0), (0, 0, 1), 1.0)
    with pytest.raises(nf.GeometryError):
        nf.hatch([circle], (1, 0), 0.0)
    with pytest.raises(nf.GeometryError):
        nf.hatch([circle], (1, 0), 0.1, contours=[1.5])
    with pytest.raises(ValueError):
        nf.hatch([circle], (1,), 0.1)
