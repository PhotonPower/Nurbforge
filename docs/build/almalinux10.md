# Building on AlmaLinux 10

These steps mirror the `almalinux10` CI job (`.github/workflows/ci.yml`, container
`almalinux:10`), which runs them on every push.

## System packages

```bash
sudo dnf install -y gcc-c++ git python3 python3-devel python3-pip
```

GCC 14 (the AlmaLinux 10 default) is used. The gmsh pip wheel needs a few X/GL runtime libraries
even when used headless:

```bash
sudo dnf install -y mesa-libGL mesa-libGLU libglvnd-opengl libgomp fontconfig \
                    libX11 libXext libXfixes libXcursor libXinerama libXft libXrender
```

Without them `import gmsh` fails with an `OSError` (missing `libGLU.so.1`, `libXft.so.2`, ...);
the mesh tests are then skipped locally, but CI treats it as an error.

## Python environment

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install --upgrade pip cmake ninja pybind11 scikit-build-core
```

CMake and Ninja come from pip, so the distribution versions do not matter.

## C++ core and tests

```bash
cmake --preset release -DNF_BUILD_PYTHON=OFF
cmake --build --preset release -j
ctest --preset release
```

Presets: `debug`, `release`, and `ci` (Release with `-Werror`). Build trees go to
`build/<preset>`.

## Python package and tests

```bash
pip install -e ".[dev]"
python -c "import gmsh; print(gmsh.__version__)"
pytest
ruff check python tests examples
```

## Wheels

Wheels built on AlmaLinux 10 target `manylinux_2_28`-compatible systems once cibuildwheel is set
up (see [roadmap](../roadmap.md)); until then, build from source as above.
