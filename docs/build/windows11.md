# Building on Windows 11

## Prerequisites

| Tool | Version | Notes |
|---|---|---|
| Visual Studio | 2022 or 2026, workload *Desktop development with C++* | Community edition is sufficient |
| CMake | ≥ 3.24; **≥ 4.2 for Visual Studio 2026** | older CMake does not know VS 2026 and falls back to NMake (error `Running 'nmake' '-?' failed`) |
| Python | ≥ 3.10 | from python.org, `winget install Python.Python.3.12`, or [uv](https://docs.astral.sh/uv/) |
| Git | any | |

The simplest way to get a matching CMake is to install it into the virtual environment (below);
it then takes precedence over an older system CMake.

## Python environment

PowerShell, from the repository root:

```powershell
python -m venv .venv                     # or: uv venv --python 3.12 .venv
.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip cmake ninja pybind11 scikit-build-core
```

With uv instead of pip: `uv pip install --python .venv cmake ninja pybind11 scikit-build-core`.

## C++ core and tests

The `msvc` preset uses the newest installed Visual Studio (multi-config generator):

```powershell
cmake --preset msvc -DNF_BUILD_PYTHON=OFF
cmake --build --preset msvc              # Release
ctest --preset msvc
```

No *Developer PowerShell* is needed – CMake locates MSVC itself. To use the single-config Ninja
presets (`debug`, `release`) instead, run them from a *Developer PowerShell for VS* so that
`cl.exe` is on `PATH`.

## Python package and tests

```powershell
pip install -e ".[dev]"                  # builds nurbforge._core with scikit-build-core
python -c "import gmsh; print(gmsh.__version__)"
pytest
ruff check python tests examples
```

`pyproject.toml` requires CMake ≥ 4.2 on Windows; if the CMake on `PATH` is older,
scikit-build-core downloads a suitable one automatically, so `pip install .` also works without
the venv CMake.

## Troubleshooting

- **`Python wurde nicht gefunden` / `Python was not found`** – the Microsoft Store alias is active
  but no Python is installed. Install Python (see above) or disable the alias under
  *Settings → Apps → Advanced app settings → App execution aliases*.
- **`Building for: NMake Makefiles`** – CMake is too old for the installed Visual Studio. Use the
  venv CMake (`pip install cmake`) or upgrade the system CMake.
- **Stale build after switching compilers/CMake** – delete `build/msvc` and configure again.
