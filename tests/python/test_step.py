"""Native STEP writer: string checks and round trips through gmsh / OpenCASCADE."""
import math

import nurbforge as nf
import pytest

try:  # gmsh raises OSError (not ImportError) if system libraries are missing
    import gmsh
except (ImportError, OSError):  # pragma: no cover
    gmsh = None

needs_gmsh = pytest.mark.skipif(gmsh is None, reason="gmsh not usable")


def sphere(c=(0, 0, 0), r=1.0):
    arc = nf.make_arc(c, (0, 0, -1), (1, 0, 0), r, 0.0, math.pi)
    return nf.Shell.sew([nf.TrimmedSurface(nf.revolve(arc, c, (0, 0, 1), 2 * math.pi))])


def box(h=1.0):  # cube [-h, h]^3
    s, rect = 2 * h, nf.make_rectangle
    return nf.Shell.sew([nf.TrimmedSurface(f) for f in [
        rect((-h, -h, -h), (s, 0, 0), (0, s, 0)), rect((-h, -h, h), (s, 0, 0), (0, s, 0)),
        rect((-h, -h, -h), (s, 0, 0), (0, 0, s)), rect((-h, h, -h), (s, 0, 0), (0, 0, s)),
        rect((-h, -h, -h), (0, s, 0), (0, 0, s)), rect((h, -h, -h), (0, s, 0), (0, 0, s))]])


def cylinder(r, z0, z1):
    side = nf.extrude(nf.make_circle((0, 0, z0), (0, 0, 1), r), (0, 0, z1 - z0))
    lids = [nf.TrimmedSurface(nf.make_rectangle((-r, -r, z), (2 * r, 0, 0), (0, 2 * r, 0)),
                              nf.circle_loop(0.5, 0.5, 0.5)) for z in (z0, z1)]
    return nf.Shell.sew([nf.TrimmedSurface(side), *lids])


def test_to_step_string():
    text = nf.to_step(box(), name="cube", timestamp="2026-01-01T00:00:00")
    assert text.startswith("ISO-10303-21;\nHEADER;\n")
    assert text.endswith("END-ISO-10303-21;\n")
    assert "AUTOMOTIVE_DESIGN" in text and "SI_UNIT(.MILLI.,.METRE.)" in text
    assert text.count("ADVANCED_FACE(") == 6
    assert text.count("EDGE_CURVE(") == 12
    assert text.count("VERTEX_POINT(") == 8
    assert nf.to_step([box()], name="cube", timestamp="2026-01-01T00:00:00") == text
    assert "CONVERSION_BASED_UNIT('INCH'" in nf.to_step(box(), unit="inch")
    with pytest.raises(ValueError):
        nf.to_step(box(), unit="furlong")


def test_write_step_file(tmp_path):
    path = tmp_path / "s.step"
    nf.write_step(path, [sphere()])
    text = path.read_text()
    assert "RATIONAL_B_SPLINE_SURFACE" in text and "MANIFOLD_SOLID_BREP" in text


def _roundtrip(tmp_path, shell, name):
    """Write natively, re-import with OCC; returns (number of volumes, OCC volume)."""
    path = str(tmp_path / f"{name}.step")
    nf.write_step(path, shell)
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        ents = gmsh.model.occ.importShapes(path)
        gmsh.model.occ.synchronize()
        vols = [t for d, t in ents if d == 3]
        return len(vols), (gmsh.model.occ.getMass(3, vols[0]) if vols else 0.0)
    except OSError as exc:  # pragma: no cover - missing runtime libraries
        pytest.skip(f"gmsh runtime not available: {exc}")
    finally:
        gmsh.finalize()


@needs_gmsh
@pytest.mark.parametrize("name", ["cube", "sphere", "union", "drilled", "void"])
def test_step_roundtrip_occ(tmp_path, name):
    shell, rel = {
        "cube": lambda: (box(), 1e-9),  # planar faces: OCC's volume is exact
        # curved faces: OCC's mass integration is approximate (observed <= 4e-4)
        "sphere": lambda: (sphere(), 2e-3),
        "union": lambda: (nf.boolean(sphere(), sphere((1, 0, 0)), nf.BooleanOp.Union), 2e-3),
        "drilled": lambda: (nf.boolean(box(), cylinder(0.3, -2, 2), nf.BooleanOp.Difference), 2e-3),
        "void": lambda: (nf.boolean(box(), sphere((0.1, -0.1, 0.2), 0.5), nf.BooleanOp.Difference),
                         2e-3),
    }[name]()
    ref = shell.oriented().volume() if not shell.is_consistently_oriented else shell.volume()
    n, vol = _roundtrip(tmp_path, shell, name)
    assert n == 1
    assert abs(vol - ref) <= rel * ref, (vol, ref)


@needs_gmsh
def test_step_roundtrip_units(tmp_path):
    # OCC converts to millimetres on import: a 2 m cube has 8e9 mm^3
    path = str(tmp_path / "m.step")
    nf.write_step(path, box(), unit="m")
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        ents = gmsh.model.occ.importShapes(path)
        gmsh.model.occ.synchronize()
        assert abs(gmsh.model.occ.getMass(3, ents[0][1]) - 8e9) < 1e-3
    finally:
        gmsh.finalize()
