"""Build the user manual PDF.

    python docs/manual/build.py [--out DIR]

1. generates the command reference from the installed ``nurbforge`` package
   (``gen_reference.py --check``: fails if an API entry has no description),
2. runs latexmk (or pdflatex three times) on ``nurbforge-manual.tex``.

The PDF is written to ``docs/manual/build/nurbforge-manual.pdf`` (or ``--out``).
Requirements: the nurbforge package, pdfLaTeX with the packages listed in ``preamble.tex``
(TeX Live: collections latex-recommended, latex-extra, fonts-recommended, fonts-extra).
"""
from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", type=Path, default=HERE / "build")
    args = ap.parse_args()
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=True)

    subprocess.run([sys.executable, str(HERE / "gen_reference.py"), "--check"], check=True)
    tex = "nurbforge-manual.tex"
    if shutil.which("latexmk"):
        cmd = ["latexmk", "-pdf", "-interaction=nonstopmode", "-halt-on-error",
               f"-outdir={out}", tex]
        subprocess.run(cmd, cwd=HERE, check=True)
    else:
        for _ in range(3):  # table of contents and cross references need three passes
            subprocess.run(["pdflatex", "-interaction=nonstopmode", "-halt-on-error",
                            f"-output-directory={out}", tex], cwd=HERE, check=True,
                           stdout=subprocess.DEVNULL)
    pdf = out / "nurbforge-manual.pdf"
    print(f"wrote {pdf}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
