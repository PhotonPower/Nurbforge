"""Check a Nurbforge installation: version, compiled core and optional components."""
import platform
import sys

import nurbforge as nf

print(f"Nurbforge {nf.__version__} on Python {platform.python_version()} "
      f"({platform.system()} {platform.machine()})")
print("compiled core:", nf._core.__file__.rsplit("/", 1)[-1].rsplit("\\", 1)[-1])

# a quick self-test of the kernel: the volume of a unit sphere
v = nf.make_sphere((0, 0, 0), 1.0).volume()
print(f"unit sphere volume {v:.12f}")

try:  # gmsh is optional (pip install "nurbforge[mesh]")
    import gmsh
    print("gmsh", gmsh.__version__)
except (ImportError, OSError) as exc:  # OSError: gmsh found, but system libraries missing
    print("gmsh not usable:", exc)

sys.exit(0 if abs(v - 4.18879020478639) < 1e-8 else 1)
