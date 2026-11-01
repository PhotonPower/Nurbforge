"""Run every Python example so the examples cannot rot.

Each script runs in a fresh interpreter inside a temporary directory (examples write mesh/STEP
files to the working directory). Scripts that mesh or export catch a missing gmsh themselves; the
second test hides gmsh to check that they still run without it. The Ollama example needs a
running server and is skipped.
"""
import re
import subprocess
import sys
from pathlib import Path

import pytest

EXAMPLES = Path(__file__).resolve().parents[2] / "examples" / "python"
NEEDS_SERVER = {"04_llm_ollama.py"}
SCRIPTS = sorted(p for p in EXAMPLES.glob("[0-9][0-9]_*.py") if p.name not in NEEDS_SERVER)
USES_GMSH = [p for p in SCRIPTS if "nurbforge.mesh" in p.read_text("utf-8")]

# runs a script with `import gmsh` failing, as on a machine without gmsh
WITHOUT_GMSH = """
import runpy, sys
sys.modules["gmsh"] = None  # makes `import gmsh` raise ImportError
runpy.run_path(sys.argv[1], run_name="__main__")
"""


def run(args, cwd):
    proc = subprocess.run([sys.executable, *args], cwd=cwd, capture_output=True, text=True,
                          timeout=120, check=False)
    assert proc.returncode == 0, f"{args[-1]} failed:\n{proc.stdout}\n{proc.stderr}"
    assert proc.stdout.strip(), f"{args[-1]} printed nothing"
    return proc.stdout


def test_examples_found():
    assert len(SCRIPTS) >= 9 and USES_GMSH


@pytest.mark.parametrize("script", SCRIPTS, ids=lambda p: p.name)
def test_example_runs(script, tmp_path):
    run([str(script)], tmp_path)


@pytest.mark.parametrize("script", USES_GMSH, ids=lambda p: p.name)
def test_example_runs_without_gmsh(script, tmp_path):
    out = run(["-c", WITHOUT_GMSH, str(script)], tmp_path)
    assert "skipped" in out
    assert not any(tmp_path.iterdir())  # nothing was written


def test_tutorial_snippets_match_examples():
    """Code blocks in docs/tutorial.md starting with `# from examples/python/<script>` must be
    verbatim excerpts of that script (every line occurs in it, in the same order)."""
    tutorial = (EXAMPLES.parents[1] / "docs" / "tutorial.md").read_text("utf-8")
    pattern = re.compile(r"```python\n# from examples/python/(\S+)\n(.*?)```", re.DOTALL)
    blocks = pattern.findall(tutorial)
    assert len(blocks) >= 10
    for name, code in blocks:
        script = (EXAMPLES / name).read_text("utf-8").splitlines()
        pos = 0
        for line in code.rstrip("\n").splitlines():
            if not line.strip():
                continue
            assert line in script[pos:], f"tutorial snippet line not in {name} (in order): {line!r}"
            pos = script.index(line, pos) + 1
