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
