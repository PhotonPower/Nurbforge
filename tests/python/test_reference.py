"""The compact agent reference (python/nurbforge/REFERENCE.md): shipped with the package,
up to date with the API (gen_agent_reference.py --check), and its example runs."""
import os
import re
import subprocess
import sys
from pathlib import Path

import nurbforge as nf

GEN = Path(__file__).resolve().parents[2] / "docs" / "manual" / "gen_agent_reference.py"


def test_reference_is_installed():
    text = nf.reference()
    assert text.startswith(f"# nurbforge {nf.__version__} ")
    for name in ("make_box", "boolean", "Shell.sew", "TrimmedSurface", "write_step", "bem_mesh"):
        assert name in text


def test_python_m_nurbforge_prints_reference():
    env = dict(os.environ, PYTHONIOENCODING="utf-8")
    res = subprocess.run([sys.executable, "-m", "nurbforge"], capture_output=True, text=True,
                         encoding="utf-8", env=env, timeout=60, check=True)
    assert res.stdout == nf.reference()


def test_reference_is_up_to_date():
    res = subprocess.run([sys.executable, str(GEN), "--check"], capture_output=True, text=True,
                         timeout=300, check=False)
    assert res.returncode == 0, res.stderr


def test_reference_example_runs(tmp_path):
    code = re.search(r"```python\n(.*?)```", nf.reference(), re.DOTALL).group(1)
    script = tmp_path / "example.py"
    script.write_text(code, encoding="utf-8")
    res = subprocess.run([sys.executable, str(script)], cwd=tmp_path, capture_output=True,
                         text=True, timeout=300, check=False)
    assert res.returncode == 0, res.stderr
    assert (tmp_path / "part.step").exists()
