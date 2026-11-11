# Building on macOS

These steps mirror the `macos` CI job (`.github/workflows/ci.yml`, `macos-latest`, Apple
Silicon), which runs them on every push. Intel Macs work the same way.

## Requirements

- macOS 14 (Sonoma) or newer
- Xcode Command Line Tools (Apple Clang with C++20): `xcode-select --install`
- Python ≥ 3.10, e.g. from [python.org](https://www.python.org/downloads/macos/) or Homebrew
  (`brew install python`)

No further system packages are needed; the gmsh pip wheel runs headless on macOS as is.

## Python environment

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install --upgrade pip cmake ninja pybind11 scikit-build-core
```

CMake and Ninja come from pip, so no Homebrew versions are needed.

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

## Deployment target

Without further settings Apple Clang builds for the running macOS version, so binaries may not
start on older systems. To build for other Macs, set the minimum version before configuring:

```bash
export MACOSX_DEPLOYMENT_TARGET=14.0
```

The release wheels use 14.0. The lower bound is 13.3, from which Apple's libc++ provides
`std::to_chars` for `double` (STEP writer); CMake stops with a clear message below that.

## Wheels

The [Wheels workflow](../../.github/workflows/wheels.yml) builds `arm64` and `x86_64` wheels on
an Apple Silicon runner (the `x86_64` wheels are tested under Rosetta 2). Locally:

```bash
pipx run cibuildwheel --only cp312-macosx_arm64
```
