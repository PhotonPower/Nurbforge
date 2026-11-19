"""Fitting surfaces to scattered points: height fields, noise, boundary curves and a base."""
import math

import numpy as np
import nurbforge as nf

rng = np.random.default_rng(1)


def sag(x, y):
    """A freeform test surface z = f(x, y)."""
    return 0.05 * (x ** 2 + 0.5 * y ** 2) + 0.02 * np.sin(2 * x) * np.cos(3 * y)


def deviation(surface, f, inner=0.05):
    """Max |z - f(x, y)| on a grid inside the domain (between and away from the data)."""
    t = np.linspace(inner, 1 - inner, 30)
    return max(abs(p.z - f(p.x, p.y)) for p in (surface(u, v) for u in t for v in t))


# --- unordered measurement of a height map over a disc of radius 2
n = 3000
r, a = 2 * np.sqrt(rng.random(n)), 2 * np.pi * rng.random(n)
x, y = r * np.cos(a), r * np.sin(a)
exact = np.column_stack([x, y, sag(x, y)])
fit = nf.fit_surface(exact)
print(f"exact data: {fit.surface.num_u} x {fit.surface.num_v} control points, "
      f"max error {fit.max_error:.1e}, rms {fit.rms_error:.1e}")
print(f"  deviation from f between the points (inner disc): "
      f"{deviation(fit.surface, sag, 0.25):.1e}")

# --- the same with measurement noise (sigma = 2 um for a part in mm)
noisy = exact + np.column_stack([np.zeros(n), np.zeros(n), 0.002 * rng.standard_normal(n)])
for num, s in ((0, 1e-3), (0, 1.0), (12, 0.1)):
    f = nf.fit_surface(noisy, num_u=num, num_v=num, smoothing=s)
    print(f"noisy data, {f.surface.num_u} x {f.surface.num_v}, smoothing {s:g}: rms residual "
          f"{f.rms_error:.2e}, deviation from f {deviation(f.surface, sag, 0.25):.1e}")

# --- boundary curves: a patch spanned by four given edges, interpolated exactly
edges = [nf.make_line((0, 0, 0), (1, 0, 0)), nf.make_line((0, 1, 0), (1, 1, 0)),
         nf.make_line((0, 0, 0), (0, 1, 0)), nf.make_line((1, 0, 0), (1, 1, 0))]
u, v = rng.random(800), rng.random(800)
bump = np.column_stack([u, v, 0.2 * np.sin(math.pi * u) * np.sin(math.pi * v)])
patch = nf.fit_surface(bump, boundary=edges)
rim = max(abs(patch.surface(t, e).z) for t in np.linspace(0, 1, 11) for e in (0.0, 1.0))
print(f"boundary curves: max error {patch.max_error:.1e}, rim height {rim:.1e}")

# --- a base surface for data that are not a height field: a bumpy half cylinder
base = nf.revolve(nf.make_line((0, 0, -1), (3, 0, -1)), (0, 0, 0), (1, 0, 0), math.pi)
x, a = 3 * rng.random(1500), math.pi * rng.random(1500)
rr = 1 + 0.05 * np.sin(2 * x) * np.sin(2 * a)
tube = np.column_stack([x, rr * np.sin(a), -rr * np.cos(a)])
fit = nf.fit_surface(tube, base=base)
print(f"half cylinder with base: max error {fit.max_error:.1e}")
