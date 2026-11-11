"""Face labels through Booleans and as gmsh physical groups."""

import math

import nurbforge as nf
import pytest


def test_labels_survive_a_boolean():
    plate = nf.make_box((-2, -2, 0), (4, 4, 1))
    tool = nf.make_cylinder((0, 0, -1), (0, 0, 1), 1.0, 3.0).with_face_labels(["bore", "", ""])
    part = nf.boolean(plate, tool, nf.BooleanOp.Difference)
    assert set(part.face_labels) == {"bottom", "top", "front", "back", "left", "right", "bore"}
    wall = sum(part.faces[i].area() for i in part.faces_labeled("bore"))
    assert math.isclose(wall, 2 * math.pi, rel_tol=1e-6)
    assert part.faces[0].with_label("x").label == "x"


def test_physical_groups_from_labels(tmp_path):
    try:  # gmsh raises OSError (not ImportError) if system libraries are missing
        import gmsh
    except (ImportError, OSError) as exc:  # pragma: no cover
        pytest.skip(f"gmsh not usable: {exc}")
    from nurbforge.mesh import GmshSession, add_shell

    part = nf.make_cylinder((0, 0, 0), (0, 0, 1), 1.0, 2.0).with_face_labels(
        ["wall", "inlet", "outlet"])
    with GmshSession():
        add_shell(part, physical_groups=True)
        groups = {gmsh.model.getPhysicalName(d, t): (d, t) for d, t in gmsh.model.getPhysicalGroups()}
        assert set(groups) == {"wall", "inlet", "outlet", "volume"}
        for name, area in (("inlet", math.pi), ("outlet", math.pi), ("wall", 4 * math.pi)):
            d, t = groups[name]
            ents = gmsh.model.getEntitiesForPhysicalGroup(d, t)
            got = sum(gmsh.model.occ.getMass(2, e) for e in ents)
            assert math.isclose(got, area, rel_tol=5e-3), name  # OCC mass integration
            # the inlet is the base disc at z = 0, the outlet at z = 2
            if name != "wall":
                z = gmsh.model.occ.getCenterOfMass(2, ents[0])[2]
                assert math.isclose(z, 0.0 if name == "inlet" else 2.0, abs_tol=1e-9)
