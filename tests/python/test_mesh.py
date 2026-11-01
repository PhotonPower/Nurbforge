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
