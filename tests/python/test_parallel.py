"""Parallel use from Python: geometry objects pickle (processes, deepcopy) and the expensive
operations release the GIL (threads)."""
import copy
import math
import multiprocessing as mp
import pickle
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor

import nurbforge as nf
import pytest


def _part():
    return nf.boolean(nf.make_box((-2, -1, 0), (4, 2, 1)),
                      nf.make_cylinder((0, 0, -1), (0, 0, 1), 0.4, 3), nf.BooleanOp.Difference)


def test_pickle_roundtrip_of_all_value_types():
    part = _part()
    objects = [
        nf.Vec3(1, 2, 3), nf.Transform.rotation((0, 0, 0), (0, 0, 1), 0.3), nf.Tolerance(1e-8),
        nf.Plane((1, 2, 3), (0, 1, 0)), nf.make_circle((0, 0, 0), (0, 0, 1), 2.0),
        nf.make_sphere((0, 0, 0), 1.0).faces[0].surface, nf.circle_loop(0.5, 0.5, 0.2),
        part.faces[0], part, nf.Asphere(curvature=0.02, conic=-1, even=[1e-6]),
        nf.ZernikeSag(base=nf.Asphere(curvature=0.01), norm_radius=5, coeffs=[0, 0, 0, 0.1]),
        nf.XYPolynomialSag(nf.Asphere(), 2.0, [(2, 0, 0.1)]), nf.SweepOptions(twist=0.3),
        nf.SectionOptions(), nf.SurfaceIntersectionOptions(chord_tolerance=1e-3),
    ]
    for o in objects:
        r = pickle.loads(pickle.dumps(o))
        assert type(r) is type(o)
    q = pickle.loads(pickle.dumps(part))
    assert q.is_closed and q.face_labels == part.face_labels
    assert q.volume() == pytest.approx(part.volume(), rel=1e-12)
    c = nf.make_circle((0, 0, 0), (0, 0, 1), 2.0)
    assert nf.distance(pickle.loads(pickle.dumps(c))(0.3), c(0.3)) == 0.0
    assert copy.deepcopy(part).volume() == pytest.approx(part.volume(), rel=1e-12)


def test_threads_give_the_same_results():
    plate = nf.make_box((-3, -2, 0), (6, 4, 1))

    def drill(r):
        return nf.boolean(plate, nf.make_cylinder((0, 0, -1), (0, 0, 1), r, 3),
                          nf.BooleanOp.Difference).volume()

    radii = [0.2 + 0.05 * i for i in range(8)]
    with ThreadPoolExecutor(4) as ex:
        par = list(ex.map(drill, radii))
    assert par == [drill(r) for r in radii]


def test_processes_return_shells():
    plate = nf.make_box((-3, -2, 0), (6, 4, 1))
    radii = [0.3, 0.5, 0.7]
    tools = [nf.make_cylinder((0, 0, -1), (0, 0, 1), r, 3) for r in radii]
    with ProcessPoolExecutor(2, mp_context=mp.get_context("spawn")) as ex:
        parts = list(ex.map(nf.boolean, [plate] * 3, tools, [nf.BooleanOp.Difference] * 3))
    for p, r in zip(parts, radii):
        assert p.is_closed
        assert p.volume() == pytest.approx(24 - math.pi * r * r, rel=1e-7)
