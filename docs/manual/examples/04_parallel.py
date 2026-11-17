"""Parameter studies in parallel: threads (the expensive operations release the GIL) and
processes (geometry objects can be pickled)."""
import math
import os
import pickle
import time
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor

import nurbforge as nf

plate = nf.make_box((-3, -2, 0), (6, 4, 1))


def drilled(r):
    """A plate with a bore of radius r - one variant of a parameter study."""
    return nf.boolean(plate, nf.make_cylinder((0, 0, -1), (0, 0, 1), r, 3), nf.BooleanOp.Difference)


if __name__ == "__main__":  # needed for process pools on macOS and Windows (spawn)
    radii = [0.2 + 0.05 * i for i in range(12)]

    t = time.perf_counter()
    serial = [drilled(r).volume() for r in radii]
    t_serial = time.perf_counter() - t

    t = time.perf_counter()
    with ThreadPoolExecutor() as ex:
        threaded = list(ex.map(lambda r: drilled(r).volume(), radii))
    t_threads = time.perf_counter() - t

    with ProcessPoolExecutor() as ex:  # the Shells come back pickled
        parts = list(ex.map(drilled, radii))

    ok = all(abs(p.volume() - (24 - math.pi * r * r)) < 1e-6 for p, r in zip(parts, radii))
    print(f"{len(radii)} variants: threads give the same volumes: {threaded == serial}, "
          f"processes return closed shells: {ok}")
    print(f"speed-up with threads on {os.cpu_count()} cores: {t_serial / t_threads:.1f}x "
          f"(machine dependent)")
    blob = pickle.dumps(parts[0])
    print(f"a pickled part: {len(blob)} bytes, faces {pickle.loads(blob).face_labels}")
