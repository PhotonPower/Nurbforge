import math
import struct

import numpy as np
import nurbforge as nf
import pytest
from nurbforge.mesh import SurfaceMesh, shrink_wrap


def _assert_closed(m, components=1, genus=0):
    assert isinstance(m, SurfaceMesh)
    assert m.order == 1
    assert m.is_watertight()
    assert m.is_consistently_oriented()
    assert np.all(m.areas() > 0)
    # Euler characteristic of closed orientable surfaces: V - E + F = 2 (components - genus)
    assert len(m.nodes) - len(m.triangles) // 2 == 2 * (components - genus)


@pytest.mark.parametrize("h, rel", [(0.2, 0.05), (0.1, 0.015)])
def test_sphere_lies_on_the_surface(h, rel):
    m = shrink_wrap(nf.make_sphere((0, 0, 0), 1.0), h)
    _assert_closed(m)
    assert m.volume() == pytest.approx(4.0 / 3.0 * math.pi, rel=rel)
    r = np.linalg.norm(m.nodes, axis=1)
    assert r.max() < 1.0 + 0.05 * h and r.min() > 1.0 - 0.15 * h
    # normals point outwards
    assert np.all(np.einsum("ij,ij->i", m.normals(), m.centroids()) > 0)


def test_offset_envelope():
    m = shrink_wrap(nf.make_sphere((0, 0, 0), 1.0), 0.1, offset=0.5)
    _assert_closed(m)
    r = np.linalg.norm(m.nodes, axis=1)
    assert r.max() == pytest.approx(1.5, abs=0.002) and r.min() == pytest.approx(1.5, abs=0.01)


def test_box_keeps_its_extent():
    m = shrink_wrap(nf.make_box((0, 0, 0), (2, 1, 1)), 0.05)
    _assert_closed(m)
    assert m.volume() == pytest.approx(2.0, rel=0.01)
    assert np.allclose(m.nodes.min(axis=0), 0.0, atol=0.005)
    assert np.allclose(m.nodes.max(axis=0), (2, 1, 1), atol=0.005)


def test_overlapping_solids_become_one_hull():
    box = nf.make_box((0, 0, 0), (2, 1, 1))
    ball = nf.make_sphere((2, 0.5, 0.5), 0.8)
    union = nf.boolean(box, ball, nf.BooleanOp.Union)
    m = shrink_wrap([box, ball], 0.05)
    _assert_closed(m)
    assert m.volume() == pytest.approx(union.volume(), rel=0.01)


def test_torus_keeps_its_hole():
    m = shrink_wrap(nf.make_torus((0, 0, 0), (0, 0, 1), 2.0, 0.5), 0.1)
    _assert_closed(m, genus=1)
    assert m.volume() == pytest.approx(2 * math.pi**2 * 2.0 * 0.25, rel=0.03)


def test_gaps_below_the_resolution_are_closed():
    a, b = nf.make_box((0, 0, 0), (1, 1, 1)), nf.make_box((1.04, 0, 0), (1, 1, 1))
    closed = shrink_wrap([a, b], 0.05)
    _assert_closed(closed)
    assert closed.volume() == pytest.approx(2.04, rel=0.01)
    apart = shrink_wrap([a, b], 0.02)
    _assert_closed(apart, components=2)
    assert apart.volume() == pytest.approx(2.0, rel=0.01)


def test_open_and_loose_surfaces():
    # a box without its lid: the opening is wide, so the hull follows the inside as well
    box = nf.make_box((0, 0, 0), (1, 1, 1))
    faces = sorted(box.faces, key=lambda f: -f.surface.evaluate(0.5, 0.5).z)
    lid, rest = faces[0], faces[1:]
    assert lid.surface.evaluate(0.5, 0.5).z == pytest.approx(1.0)
    cup = shrink_wrap(rest, 0.05, offset=0.05)
    _assert_closed(cup)
    assert cup.volume() == pytest.approx(1.1**2 * 1.05 - 0.9**2 * 0.95, rel=0.03)
    # with the lid (as a loose surface) the box is closed again
    full = shrink_wrap([*rest, lid.surface], 0.05, offset=0.05)
    _assert_closed(full)
    assert full.volume() > 1.25

    sheet = nf.make_bilinear((0, 0, 0), (1, 0, 0), (0, 1, 0), (1, 1, 0))
    m = shrink_wrap(sheet, 0.05, offset=0.1)
    _assert_closed(m)
    exact = 0.2 + 4 * 0.5 * math.pi * 0.01 + 4.0 / 3.0 * math.pi * 1e-3
    assert m.volume() == pytest.approx(exact, rel=0.03)
    with pytest.raises(ValueError, match="nothing is left"):
        shrink_wrap(sheet, 0.05)


def test_fill_holes():
    ball = nf.make_sphere((0, 0, 0), 1.0)
    m = shrink_wrap(ball, 0.1, offset=0.2, fill_holes=False)
    _assert_closed(m, components=2)
    assert m.volume() == pytest.approx(4.0 / 3.0 * math.pi * (1.2**3 - 0.8**3), rel=0.03)
    assert shrink_wrap(ball, 0.1, offset=0.2).volume() == pytest.approx(
        4.0 / 3.0 * math.pi * 1.2**3, rel=0.02)


def test_trimmed_faces():
    part = nf.boolean(nf.make_box((-1, -1, -1), (2, 2, 2)),
                      nf.make_cylinder((0, 0, -2), (0, 0, 1), 0.5, 4.0),
                      nf.BooleanOp.Difference)
    m = shrink_wrap(part, 0.05)
    _assert_closed(m, genus=1)
    assert m.volume() == pytest.approx(part.volume(), rel=0.015)


def test_errors():
    box = nf.make_box((0, 0, 0), (1, 1, 1))
    with pytest.raises(ValueError, match="resolution"):
        shrink_wrap(box, 0.0)
    with pytest.raises(ValueError, match="offset"):
        shrink_wrap(box, 0.1, offset=-0.1)
    with pytest.raises(ValueError, match="max_nodes"):
        shrink_wrap(box, 0.001)
    with pytest.raises(ValueError, match="nothing to wrap"):
        shrink_wrap([], 0.1)
    with pytest.raises(TypeError):
        shrink_wrap([nf.make_line((0, 0, 0), (1, 0, 0))], 0.1)
    with pytest.raises(TypeError):
        shrink_wrap(3.0, 0.1)


def test_write(tmp_path):
    m = shrink_wrap(nf.make_box((0, 0, 0), (1, 1, 1)), 0.2)
    stl = tmp_path / "hull.stl"
    m.write(stl)
    data = stl.read_bytes()
    (count,) = struct.unpack_from("<I", data, 80)
    assert count == len(m.triangles) and len(data) == 84 + 50 * count
    first = np.frombuffer(data, dtype="<f4", count=12, offset=84).reshape(4, 3)
    assert np.allclose(first[0], m.normals()[0], atol=1e-6)
    assert np.allclose(first[1:], m.nodes[m.triangles[0]], atol=1e-6)

    msh = tmp_path / "hull.msh"
    m.write(msh)
    lines = msh.read_text().splitlines()
    assert lines[:3] == ["$MeshFormat", "2.2 0 8", "$EndMeshFormat"]
    assert int(lines[lines.index("$Nodes") + 1]) == len(m.nodes)
    e = lines.index("$Elements")
    assert int(lines[e + 1]) == len(m.triangles)
    assert [int(x) for x in lines[e + 2].split()][5:] == (m.triangles[0] + 1).tolist()
    with pytest.raises(ValueError, match="unsupported"):
        m.write(tmp_path / "hull.xyz")


def test_written_msh_is_read_by_gmsh(tmp_path):
    gmsh = pytest.importorskip("gmsh")
    m = shrink_wrap(nf.make_sphere((0, 0, 0), 1.0), 0.2)
    path = tmp_path / "hull.msh"
    m.write(path)
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.open(str(path))
        tags, _, _ = gmsh.model.mesh.getNodes()
        _, elements, _ = gmsh.model.mesh.getElements(2)
        assert len(tags) == len(m.nodes) and len(elements[0]) == len(m.triangles)
    finally:
        gmsh.finalize()
