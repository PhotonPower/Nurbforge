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
