"""Surface analysis: curvature and draft angle maps, direction, continuity, duplicates."""
import math

import numpy as np
import nurbforge as nf
import pytest


def test_curvature_map():
    torus = nf.make_torus((0, 0, 0), (0, 0, 1), 3.0, 1.0).faces[0].surface
    m = nf.curvature_map(torus, samples_per_span=4)
    assert m.gaussian.shape == (len(m.u), len(m.v)) == m.inside.shape
    assert m.inside.all()
    assert m.gaussian_range.min == pytest.approx(-0.5, abs=1e-9)
    assert m.gaussian_range.max == pytest.approx(0.25, abs=1e-9)
    assert m.gaussian.min() == pytest.approx(m.gaussian_range.min)
    assert m.min_radius == pytest.approx(1.0, abs=1e-9)
    u, v = m.gaussian_range.min_uv
    assert math.hypot(torus(u, v).x, torus(u, v).y) == pytest.approx(2.0, abs=1e-9)
    assert np.allclose(m.k_max * m.k_min, m.gaussian, atol=1e-9)
    plate = nf.boolean(nf.make_box((0, 0, 0), (4, 2, 1)),
                       nf.make_cylinder((2, 1, -1), (0, 0, 1), 0.5, 3), nf.BooleanOp.Difference)
    trimmed = [nf.curvature_map(f, 4) for f in plate.faces]
    assert any(not t.inside.all() for t in trimmed)  # a face with the hole


def test_draft_angle():
    cone = nf.extrude_to_point(nf.make_circle((0, 0, 0), (0, 0, 1), 1.0), (0, 0, 2))
    m = nf.draft_angle_map(cone, (0, 0, 1), 4)
    assert abs(m.range.min) == pytest.approx(math.atan(0.5), abs=1e-9)
    assert abs(m.range.max) == pytest.approx(math.atan(0.5), abs=1e-9)
    assert m.angle.shape == (len(m.u), len(m.v))
    assert nf.draft_angle(cone, 0.3, 0.4, (0, 0, 1)) == pytest.approx(m.range.min, abs=1e-9)
    with pytest.raises(nf.GeometryError):
        nf.draft_angle(cone, 0.3, 0.4, (0, 0, 0))


def test_surface_direction():
    plane = nf.make_rectangle((0, 0, 0), (2, 0, 0), (0, 3, 0))
    d = nf.surface_direction(plane)
    assert nf.distance(d.normal, nf.Vec3(0, 0, 1)) < 1e-12 and d.planar
    assert nf.distance(nf.surface_direction(plane.transposed()).normal, nf.Vec3(0, 0, -1)) < 1e-12
    sphere = nf.make_sphere((0, 0, 0), 1.0).faces[0].surface
    s = nf.surface_direction(sphere)
    assert sum(s.poles) == 2 and (s.closed_u != s.closed_v) and not s.planar


def test_continuity_of_surfaces_and_shell_edges():
    floor = nf.make_rectangle((-2, 0, 0), (2, 0, 0), (0, 3, 0))
    bend = nf.extrude(nf.make_arc((0, 0, 1), (0, 0, -1), (1, 0, 0), 1.0, 0.0, math.pi / 2),
                      (0, 3, 0))
    k = nf.continuity(floor, nf.SurfaceSide.UMax, bend)
    assert k.geometric == 1 and k.max_gap < 1e-12
    assert k.max_curvature_difference == pytest.approx(1.0, abs=1e-9)
    wall = nf.make_rectangle((0, 0, 0), (0, 0, 2), (0, 3, 0))
    assert nf.continuity(floor, nf.SurfaceSide.UMax, wall).geometric == 0
    assert nf.continuity(floor, nf.SurfaceSide.UMin, wall).geometric == -1
    cyl = nf.make_cylinder((0, 0, 0), (0, 0, 1), 1.0, 2.0)
    assert sorted(e.geometric for e in nf.edge_continuity(cyl)) == [0, 0, 2]
    assert all(e.geometric == 0 for e in nf.edge_continuity(nf.make_box((0, 0, 0), (1, 2, 3))))


def test_duplicate_surfaces():
    plane = nf.make_rectangle((0, 0, 0), (2, 0, 0), (0, 3, 0))
    half = nf.make_rectangle((0, 0, 0), (1, 0, 0), (0, 3, 0))
    assert nf.find_duplicate_surfaces([plane, half, plane.transposed()], 1e-9) == [[0, 2]]
    assert nf.find_duplicate_surfaces([plane, half]) == []
