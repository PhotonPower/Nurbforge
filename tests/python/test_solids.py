"""Solid primitives, mass properties and transformed shells from Python."""

import math

import numpy as np
import nurbforge as nf
import pytest


def test_primitives_volumes():
    cases = [
        (nf.make_box((1, 2, 3), (2, 3, 4)), 24.0),
        (nf.make_sphere((0, 0, 0), 2.0), 4 / 3 * math.pi * 8),
        (nf.make_cylinder((0, 0, 0), (1, 0, 0), 1.0, 3.0), 3 * math.pi),
        (nf.make_cone((0, 0, 0), (0, 0, 1), 2.0, 0.0, 3.0), 4 * math.pi),
        (nf.make_torus((0, 0, 0), (0, 0, 1), 3.0, 1.0), 6 * math.pi**2),
    ]
    for shell, volume in cases:
        assert shell.is_closed
        assert math.isclose(shell.volume(), volume, rel_tol=1e-9)


def test_mass_properties_of_a_cylinder():
    r, h = 1.5, 4.0
    mp = nf.make_cylinder((0, 0, -1), (0, 0, 1), r, h).mass_properties()
    v = math.pi * r * r * h
    assert math.isclose(mp.volume, v, rel_tol=1e-9)
    assert np.allclose(tuple(mp.centroid), (0, 0, 1), atol=1e-9)
    expected = np.diag([v * (3 * r * r + h * h) / 12] * 2 + [0.5 * v * r * r])
    assert np.allclose(mp.inertia, expected, rtol=1e-8, atol=1e-8)
    assert math.isclose(mp.principal_moments[0], 0.5 * v * r * r, rel_tol=1e-8)
    assert np.allclose(np.abs(mp.principal_axes[0]), (0, 0, 1), atol=1e-8)


def test_transformed_shells():
    box = nf.make_box((0, 0, 0), (1, 2, 3))
    t = nf.Transform.rotation((0, 0, 0), (0, 0, 1), math.pi / 2)
    moved = t * box
    assert math.isclose(moved.volume(), 6.0, rel_tol=1e-10)
    assert np.allclose(tuple(moved.mass_properties().centroid), (-1.0, 0.5, 1.5), atol=1e-10)
    mirrored = box.transformed(nf.Transform.mirror((0, 0, 0), (0, 0, 1)))
    assert math.isclose(mirrored.volume(), 6.0, rel_tol=1e-10)  # still outward


def test_errors():
    with pytest.raises(nf.GeometryError):
        nf.make_box((0, 0, 0), (1, -1, 1))
    with pytest.raises(nf.GeometryError):
        nf.make_torus((0, 0, 0), (0, 0, 1), 1.0, 2.0)
    with pytest.raises(nf.GeometryError):
        nf.Shell.sew([nf.TrimmedSurface(nf.make_rectangle((0, 0, 0), (1, 0, 0), (0, 1, 0)))]
                     ).mass_properties()
