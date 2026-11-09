"""Curve constructors: circles through points, polygons, helices, rounded corners."""

import math

import numpy as np
import nurbforge as nf


def test_circle_and_arc_through_points():
    c = nf.make_circle_3pt((1, 0, 0), (0, 1, 0), (-1, 0, 0))
    assert np.allclose(np.linalg.norm(c.sample(40), axis=1), 1.0, atol=1e-12)
    arc = nf.make_arc_3pt((1, 0, 0), (0, 1, 0), (-1, 0, 0))
    assert math.isclose(nf.length(arc), math.pi, rel_tol=1e-10)


def test_polygon_rectangle_and_round_corners():
    square = nf.make_rectangle_curve((0, 0, 0), (2, 0, 0), (0, 2, 0))
    rounded = nf.fillet_corners(square, 0.3)
    assert math.isclose(nf.length(rounded), 8 - 8 * 0.3 + 2 * math.pi * 0.3, rel_tol=1e-10)
    tri = nf.make_polygon((0, 0, 0), (0, 0, 1), 1.0, 3)
    assert math.isclose(nf.length(tri), 3 * math.sqrt(3), rel_tol=1e-12)


def test_helix_follows_the_exact_curve():
    h = nf.make_helix((0, 0, 0), (0, 0, 1), (2, 0, 0), 0.5, 3.0, tolerance=1e-9)
    t = np.linspace(0, 3, 301)
    exact = np.stack([2 * np.cos(2 * np.pi * t), 2 * np.sin(2 * np.pi * t), 0.5 * t], axis=1)
    got = np.array([tuple(h(x)) for x in t])
    assert np.max(np.linalg.norm(got - exact, axis=1)) < 1e-9
