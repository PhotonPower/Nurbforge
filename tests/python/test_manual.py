"""The user manual (docs/manual): every example script runs, and the generated command reference
covers the whole public API (gen_reference.py --check)."""
import os
import subprocess
import sys
from pathlib import Path

import pytest

MANUAL = Path(__file__).resolve().parents[2] / "docs" / "manual"
EXAMPLES = sorted((MANUAL / "examples").glob("*.py"))


def _usable(module: str) -> bool:
    try:
        __import__(module)
    except (ImportError, OSError):  # gmsh raises OSError if its system libraries are missing
        return False
    return True


@pytest.mark.parametrize("script", EXAMPLES, ids=[p.name for p in EXAMPLES])
def test_manual_example_runs(script, tmp_path):
    first = script.read_text(encoding="utf-8").splitlines()[0]
    if first.startswith("# requires:"):
        need = first.split(":", 1)[1].strip()
        if need == "ollama":
            pytest.skip("needs a running Ollama server")
        if need == "gmsh" and not _usable("gmsh"):
            pytest.skip("gmsh not usable")
    env = dict(os.environ, PYTHONIOENCODING="utf-8", MPLBACKEND="Agg")
    res = subprocess.run([sys.executable, str(script)], cwd=tmp_path, env=env,
                         capture_output=True, text=True, timeout=300, check=False)
    assert res.returncode == 0, f"{script.name} failed:\n{res.stdout}\n{res.stderr}"


def test_command_reference_is_complete(tmp_path):
    res = subprocess.run([sys.executable, str(MANUAL / "gen_reference.py"), "--check",
                          "--out", str(tmp_path)], capture_output=True, text=True, timeout=300, check=False)
    assert res.returncode == 0, res.stderr
    ref = (tmp_path / "reference.tex").read_text(encoding="utf-8")
    assert "\\begin{refentry}{make_box}" in ref
    assert "\\begin{refentry}{mesh.bem_mesh}" in ref
    assert (tmp_path / "llm-tools.tex").read_text(encoding="utf-8").count("\\begin{toolentry}") > 30
