import math

import numpy as np
import nurbforge as nf
import pytest

try:  # gmsh raises OSError (not ImportError) if system libraries are missing
    import gmsh
except (ImportError, OSError) as exc:  # pragma: no cover
    pytest.skip(f"gmsh not usable: {exc}", allow_module_level=True)

from nurbforge.mesh import GmshSession


def test_surface_transfer_is_exact():
    arc = nf.make_arc((0, 0, 0), (1, 0, 0), (0, 1, 0), 2.0, 0.0, 1.2)
    s = nf.extrude(arc, (0.3, 0.1, 4.0))
    with GmshSession() as g:
        tag = g.add_surface(s)
        gmsh.model.occ.synchronize()
        for u, v in [(0.2, 0.3), (0.7, 0.9)]:
            assert np.allclose(gmsh.model.getValue(2, tag, [u, v]), tuple(s(u, v)), atol=1e-9)


def test_curve_transfer_is_exact():
    c = nf.make_circle((1, 1, 1), (1, 1, 0), 3.0)
    with GmshSession() as g:
        tag = g.add_curve(c)
        gmsh.model.occ.synchronize()
        assert np.allclose(gmsh.model.getValue(1, tag, [0.37]), tuple(c(0.37)), atol=1e-9)


def test_sphere_mesh_nodes_on_surface(tmp_path):
    arc = nf.make_arc((0, 0, 0), (0, 0, -1), (1, 0, 0), 1.0, 0.0, math.pi)
    sph = nf.revolve(arc, (0, 0, 0), (0, 0, 1), 2 * math.pi)
    with GmshSession() as g:
        g.add_surface(sph)
        g.mesh(dim=2, size=0.2)
        nodes = gmsh.model.mesh.getNodes()[1].reshape(-1, 3)
        g.write(str(tmp_path / "sphere.msh"))
    assert np.allclose(np.linalg.norm(nodes, axis=1), 1.0, atol=1e-9)
    assert (tmp_path / "sphere.msh").exists()


def test_trimmed_surface_mesh(tmp_path):
    plate = nf.make_rectangle((0, 0, 0), (4, 0, 0), (0, 2, 0))
    disk = nf.TrimmedSurface(plate, nf.circle_loop(0.5, 0.5, 0.25))  # ellipse, area pi/2
    try:
        with GmshSession() as g:
            tag = g.add_surface(disk)
            g.gmsh.model.occ.synchronize()
            assert abs(g.gmsh.model.occ.getMass(2, tag) - math.pi / 2) < 1e-2
            g.mesh(dim=2, size=0.1)
            g.write(str(tmp_path / "disk.msh"))
            # every mesh node lies inside the ellipse (x-2)^2 + 4 (y-1)^2 <= 1, up to the
            # tolerance of OCC's own 3D approximation of the trimmed edge (~1e-4)
            _, coords, _ = gmsh.model.mesh.getNodes()
            for i in range(0, len(coords), 3):
                x, y = coords[i], coords[i + 1]
                assert (x - 2) ** 2 + 4 * (y - 1) ** 2 <= 1 + 1e-3
    except OSError as exc:  # pragma: no cover - missing GL runtime libraries
        pytest.skip(f"gmsh runtime not available: {exc}")


def test_shell_step_export_roundtrip(tmp_path):
    R, H = 1.0, 2.0
    side = nf.extrude(nf.make_circle((0, 0, 0), (0, 0, 1), R), (0, 0, H))
    lid = [nf.TrimmedSurface(nf.make_rectangle((-R, -R, z), (2 * R, 0, 0), (0, 2 * R, 0)),
                             nf.circle_loop(0.5, 0.5, 0.5)) for z in (0.0, H)]
    shell = nf.Shell.sew([nf.TrimmedSurface(side), *lid]).oriented()
    assert shell.is_closed
    assert abs(shell.volume() - math.pi * R * R * H) < 1e-7
    path = str(tmp_path / "cylinder.step")
    try:
        from nurbforge.mesh import export_step

        export_step([shell], path)
        with GmshSession() as g:
            ents = g.gmsh.model.occ.importShapes(path)
            g.gmsh.model.occ.synchronize()
            vols = [t for d, t in ents if d == 3]
            assert len(vols) == 1
            assert len(g.gmsh.model.getEntities(2)) == 3
            # volume of a tetrahedral mesh of the re-imported solid (OCC's getMass is only
            # accurate to ~0.5% for these faces, the mesh is much closer)
            g.gmsh.option.setNumber("Mesh.MeshSizeMax", 0.1)
            g.gmsh.model.mesh.generate(3)
            _, _, conn = g.gmsh.model.mesh.getElements(3)
            ids, xyz, _ = g.gmsh.model.mesh.getNodes()
            pos = dict(zip(ids, np.asarray(xyz).reshape(-1, 3)))
            vol = 0.0
            for a, b, c, d in np.asarray(conn[0]).reshape(-1, 4):
                vol += abs(np.dot(pos[b] - pos[a], np.cross(pos[c] - pos[a], pos[d] - pos[a]))) / 6
            assert abs(vol - math.pi * R * R * H) < 3e-3 * math.pi * R * R * H
    except OSError as exc:  # pragma: no cover - missing GL runtime libraries
        pytest.skip(f"gmsh runtime not available: {exc}")


# ------------------------------------------------------------------ local mesh sizes


def _edge_lengths(tags, where=None):
    """Lengths of the triangle edges on the gmsh surfaces ``tags`` (optionally only edges whose
    midpoint satisfies ``where``)."""
    ids, xyz, _ = gmsh.model.mesh.getNodes()
    pos = dict(zip(ids, np.asarray(xyz).reshape(-1, 3)))
    out = []
    for t in tags:
        _, _, conn = gmsh.model.mesh.getElements(2, t)
        for tri in np.asarray(conn[0]).reshape(-1, 3):
            for a, b in ((0, 1), (1, 2), (2, 0)):
                p, q = pos[tri[a]], pos[tri[b]]
                if where is None or where(0.5 * (p + q)):
                    out.append(np.linalg.norm(p - q))
    return np.asarray(out)


def _labelled_box():
    return nf.make_box((0, 0, 0), (4, 4, 4))  # faces bottom/top/front/back/left/right


def test_refine_by_label():
    with GmshSession() as g:
        g.add_shell(_labelled_box())
        g.refine("top", 0.1)
        g.mesh(dim=3, size=1.0)
        tags = g.face_tags()
        top, bottom = _edge_lengths(tags["top"]), _edge_lengths(tags["bottom"])
    assert set(tags) == {"bottom", "top", "front", "back", "left", "right"}
    assert np.median(top) == pytest.approx(0.1, rel=0.15)
    assert np.median(bottom) > 0.4  # far from the refined face: global size
    assert np.max(top) < 0.2


def test_refine_with_distance_grows_towards_global_size():
    with GmshSession() as g:
        g.add_shell(_labelled_box())
        g.refine(["top"], 0.1, distance=2.0)
        g.mesh(dim=2, size=1.0)
        side = g.face_tags()["front"]
        near = _edge_lengths(side, lambda m: m[2] > 3.8)   # within 0.2 of the top
        mid = _edge_lengths(side, lambda m: 2.4 < m[2] < 2.8)
        far = _edge_lengths(side, lambda m: m[2] < 1.0)
    # linear growth: ~0.1 at the top, ~0.1 + 0.9 * 1.4 / 2 ~ 0.7 at z = 2.6, 1.0 beyond
    assert np.median(near) < 0.25
    assert 0.35 < np.median(mid) < 0.85
    assert np.median(far) > np.median(mid)


def test_refine_by_surface_tag_and_several_refinements():
    with GmshSession() as g:
        g.add_shell(_labelled_box())
        tags = g.face_tags()
        g.refine(tags["left"][0], 0.15)          # gmsh surface tag
        g.refine(["right", "top"], 0.2)          # several labels in one refinement
        g.mesh(dim=2, size=1.0)
        tags = g.face_tags()
        sizes = {k: np.median(_edge_lengths(v)) for k, v in tags.items()}
        # bottom borders the refined left/right faces, whose fine edges grade into it; its middle
        # (x, y in 1.5..2.5) stays coarse
        middle = _edge_lengths(tags["bottom"], lambda m: 1.5 < m[0] < 2.5 and 1.5 < m[1] < 2.5)
    assert sizes["left"] == pytest.approx(0.15, rel=0.15)
    assert sizes["right"] == pytest.approx(0.2, rel=0.15)
    assert sizes["top"] == pytest.approx(0.2, rel=0.15)
    assert np.median(middle) > 0.3


def test_face_tags_survive_renumbering_and_remeshing():
    # two shells: the second add_shell and mesh(heal=True) run removeAllDuplicates again
    a = _labelled_box().with_face_labels(["a"] * 6)
    b = nf.make_sphere((10, 0, 0), 1.0).with_face_labels(["ball"])
    with GmshSession() as g:
        g.add_shell(a)
        g.add_shell(b)
        g.refine("ball", 0.1)
        g.mesh(dim=2, size=1.0)
        g.mesh(dim=2, size=1.0)  # fields are rebuilt, not accumulated
        tags = g.face_tags()
        ball = _edge_lengths(tags["ball"])
        nodes = gmsh.model.mesh.getNodes(2, tags["ball"][0], includeBoundary=True)[1]
    assert len(tags["a"]) == 6 and len(tags["ball"]) == 1
    assert np.median(ball) == pytest.approx(0.1, rel=0.2)
    pts = np.asarray(nodes).reshape(-1, 3)
    assert np.allclose(np.linalg.norm(pts - (10, 0, 0), axis=1), 1.0, atol=1e-6)


def test_curvature_based_size():
    sph = nf.make_sphere((0, 0, 0), 1.0)
    with GmshSession() as g:
        g.add_shell(sph)
        g.mesh(dim=2, size=10.0, curvature=24)
        lengths = _edge_lengths(g.face_tags()["surface"])
    # 24 elements per 2*pi on a unit sphere: edges ~ 2*pi/24 ~ 0.26
    assert np.median(lengths) == pytest.approx(2 * math.pi / 24, rel=0.3)


def test_refine_errors():
    with GmshSession() as g:
        g.add_shell(_labelled_box())
        with pytest.raises(ValueError):
            g.refine("top", 0.0)
        with pytest.raises(ValueError):
            g.refine("top", 0.1, distance=-1.0)
        with pytest.raises(ValueError):
            g.refine([], 0.1)
        with pytest.raises(ValueError):
            g.mesh(dim=2, curvature=0)
        g.refine("lid", 0.1)
        with pytest.raises(KeyError, match="top"):  # message lists the existing labels
            g.mesh(dim=2, size=1.0)
