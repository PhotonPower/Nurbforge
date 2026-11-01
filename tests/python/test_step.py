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


# ---------------------------------------------------------------------------- reading


def _vol(shell):
    return abs(shell.volume(1e-7) if shell.is_consistently_oriented else shell.oriented().volume(1e-7))


@pytest.mark.parametrize("name", ["cube", "sphere", "cylinder", "void", "open"])
def test_read_step_roundtrip(tmp_path, name):
    def open_box():
        return nf.Shell.sew(list(box().faces)[:5])

    shell = {"cube": box, "sphere": lambda: sphere((0.5, 0, 1), 2.0), "cylinder": lambda: cylinder(0.5, 0, 2),
             "void": lambda: nf.boolean(box(), sphere((0.1, -0.1, 0.2), 0.5), nf.BooleanOp.Difference),
             "open": open_box}[name]()
    path = tmp_path / f"{name}.step"
    nf.write_step(path, shell)
    shells = nf.read_step(path)
    assert len(shells) == 1
    s = shells[0]
    assert len(s.faces) == len(shell.faces)
    assert s.is_closed == shell.is_closed
    assert s.area(1e-7) == pytest.approx(shell.area(1e-7), rel=1e-6)
    if shell.is_closed:
        assert s.volume(1e-7) == pytest.approx(_vol(shell), rel=1e-6)


def test_from_step_units_and_errors():
    text = nf.to_step(box(), unit="m")
    (s,) = nf.from_step(text)
    assert s.volume(1e-7) == pytest.approx(8e9, rel=1e-9)
    (s,) = nf.from_step(text, unit="m")
    assert s.volume(1e-7) == pytest.approx(8.0, rel=1e-9)
    assert nf.from_step(nf.to_step(box(), unit="inch"), unit="inch")[0].volume(1e-7) == pytest.approx(8.0)
    with pytest.raises(ValueError):
        nf.from_step(text, unit="furlong")
    with pytest.raises(nf.GeometryError):
        nf.from_step("ISO-10303-21;\nDATA;\n#1=FOO(;\n")
    bad = text.replace("B_SPLINE_SURFACE_WITH_KNOTS(''", "OFFSET_SURFACE(''", 1)
    with pytest.raises(nf.GeometryError, match="OFFSET_SURFACE"):
        nf.from_step(bad)
    with pytest.warns(UserWarning, match="OFFSET_SURFACE"):
        (s,) = nf.from_step(bad, skip_unsupported=True)
    assert len(s.faces) == 5 and not s.is_closed


def _occ_shape(name):
    """Builds a shape with gmsh's OpenCASCADE kernel; returns its exact volume (None: unknown)."""
    occ, pi = gmsh.model.occ, math.pi
    if name == "box":
        occ.addBox(0, 0, 0, 10, 20, 30)
        return 6000.0
    if name == "cylinder":
        occ.addCylinder(0, 0, 0, 0, 0, 20, 5)
        return pi * 25 * 20
    if name == "sphere":
        occ.addSphere(1, 2, 3, 4)
        return 4 / 3 * pi * 64
    if name == "cone":  # frustum
        occ.addCone(0, 0, 0, 0, 0, 10, 10, 3)
        return pi * 10 / 3 * (100 + 30 + 9)
    if name == "torus":
        occ.addTorus(0, 0, 0, 20, 5)
        return 2 * pi * pi * 20 * 25
    if name == "drilled":  # partial cylindrical faces, holes in planar faces
        occ.cut([(3, occ.addBox(0, 0, 0, 10, 10, 10))], [(3, occ.addCylinder(5, 5, -1, 0, 0, 12, 3))])
        return 1000 - pi * 9 * 10
    if name == "fillet":  # rounded box: partial cylinders, spheres
        b = occ.addBox(0, 0, 0, 10, 10, 10)
        occ.synchronize()
        occ.fillet([b], [t for _, t in gmsh.model.getEntities(1)], [2.0])
        return None
    if name == "hyperbola":  # cone cut parallel to its axis: HYPERBOLA edges
        occ.cut([(3, occ.addCone(0, 0, 0, 0, 0, 6, 3, 1))], [(3, occ.addBox(-1, -1, 2, 5, 5, 2))])
        return None
    if name == "meridian":  # sphere cut through its poles: boundary through the poles
        occ.cut([(3, occ.addSphere(0, 0, 0, 3))], [(3, occ.addBox(-4, -4, -4, 4, 8, 8))])
        return 2 / 3 * pi * 27
    if name == "void":  # BREP_WITH_VOIDS
        occ.cut([(3, occ.addBox(0, 0, 0, 10, 10, 10))], [(3, occ.addSphere(5, 5, 5, 2))])
        return 1000 - 4 / 3 * pi * 8
    raise ValueError(name)


@needs_gmsh
@pytest.mark.parametrize("name", ["box", "cylinder", "sphere", "cone", "torus", "drilled", "fillet", "hyperbola",
                                  "meridian", "void"])
def test_read_step_occ(tmp_path, name):
    """Files written by gmsh / OpenCASCADE (analytic surfaces, SURFACE_CURVE / SEAM_CURVE edges)."""
    path = str(tmp_path / f"{name}.step")
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add(name)
        exact = _occ_shape(name)
        gmsh.model.occ.synchronize()
        vols = gmsh.model.getEntities(3)
        occ_volume = sum(gmsh.model.occ.getMass(3, t) for _, t in vols)
        gmsh.write(path)
    except OSError as exc:  # pragma: no cover - missing runtime libraries
        pytest.skip(f"gmsh runtime not available: {exc}")
    finally:
        gmsh.finalize()
    shells = nf.read_step(path)
    assert len(shells) == len(vols)
    for s in shells:
        assert s.is_closed and s.is_consistently_oriented
    vol = sum(s.volume(1e-7) for s in shells)
    # OCC integrates the analytic faces of these models accurately
    assert vol == pytest.approx(exact if exact is not None else occ_volume, rel=1e-6)
