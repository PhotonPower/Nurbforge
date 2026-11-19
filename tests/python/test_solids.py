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


def test_fillet_and_chamfer_edges():
    box = nf.make_box((0, 0, 0), (4, 3, 2))
    e = nf.edges_between(box, "top", "front")
    assert len(e) == 1
    r = 0.5
    f = nf.fillet_edges(box, e, r)
    assert f.is_closed and "fillet" in f.face_labels
    assert f.volume() == pytest.approx(24 - 4 * r * r * (1 - math.pi / 4), abs=1e-7)
    c = nf.chamfer_edges(box, e, r, label="bevel")
    assert "bevel" in c.face_labels and c.volume() == pytest.approx(23.5, abs=1e-9)
    # three edges at a corner: spherical corner patch (rounded corner instead of sharp)
    corner = e + nf.edges_between(box, "top", "left") + nf.edges_between(box, "front", "left")
    k = nf.fillet_edges(box, corner, r)
    a = r * r * (1 - math.pi / 4)
    assert k.is_closed and len(k.faces_labeled("fillet")) == 4
    assert k.volume() == pytest.approx(24 - 9 * a + r ** 3 * (2 - 7 * math.pi / 12), abs=1e-7)
    # variable radius along the edge: 0.3 -> 0.6 (smooth step)
    v = nf.fillet_edges(box, e, 1.0, sizes=[((0, 0, 2), 0.3), ((4, 0, 2), 0.6)])
    integral = 4 * (0.09 + 0.3 * 0.3 + 0.09 * 181 / 462)
    assert v.volume() == pytest.approx(24 - (1 - math.pi / 4) * integral, abs=1e-7)
    with pytest.raises(nf.GeometryError):
        nf.fillet_edges(box, e, r, sizes=[((2, 2, 2), 0.3)])  # point not on the edge
    cyl = nf.make_cylinder((0, 0, 0), (0, 0, 1), 2.0, 3.0)
    rim = nf.chamfer_edges(cyl, nf.edges_between(cyl, "top", "side"), 0.4)
    assert rim.volume() == pytest.approx(12 * math.pi - 2 * math.pi * (2 - 0.4 / 3) * 0.08, abs=1e-6)
