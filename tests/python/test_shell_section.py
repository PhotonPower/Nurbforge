"""Shell sections, critical heights and named STEP shells from Python."""

import math

import nurbforge as nf


def box(h: float) -> nf.Shell:
    s = 2 * h
    rect = nf.make_rectangle
    return nf.Shell.sew([
        nf.TrimmedSurface(rect((-h, -h, -h), (s, 0, 0), (0, s, 0))),
        nf.TrimmedSurface(rect((-h, -h, -h), (s, 0, 0), (0, 0, s))),
        nf.TrimmedSurface(rect((-h, h, -h), (s, 0, 0), (0, 0, s))),
        nf.TrimmedSurface(rect((-h, -h, -h), (0, s, 0), (0, 0, s))),
        nf.TrimmedSurface(rect((h, -h, -h), (0, s, 0), (0, 0, s))),
        nf.TrimmedSurface(rect((-h, -h, h), (s, 0, 0), (0, s, 0))),
    ])


def test_section_and_heights():
    b = box(1.0)
    sec = nf.section(b, nf.Plane((0, 0, 0.25), (0, 0, 1)))
    assert len(sec.loops) == 1
    assert math.isclose(sec.loops[0].area, 4.0, rel_tol=1e-14)
    assert all(c.exact for c in sec.loops[0].curves)
    hs = nf.critical_heights(b, (0, 0, 1))
    assert len(hs) == 2 and math.isclose(hs[0], -1.0)
    assert nf.is_extrusion_band(b, (0, 0, 1), -1.0, 1.0)


def test_named_micrometre_step():
    text = nf.to_step([box(1.0), box(3.0)], unit="um", names=["Pillar", "Cap"])
    named = nf.from_step_named(text, unit="um")
    assert [n for n, _ in named] == ["Pillar", "Cap"]
    assert math.isclose(abs(named[0][1].volume()), 8.0, rel_tol=1e-9)
