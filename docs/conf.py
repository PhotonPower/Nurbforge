"""Sphinx configuration for the Nurbforge documentation site.

The Markdown pages in ``docs/`` are the single source: they render on GitHub as they are and are
parsed by MyST here. The C++ reference comes from Doxygen XML via Breathe, the Python reference
from autodoc on the installed ``nurbforge`` package.

Build locally from the repository root::

    pip install -e ".[docs]"          # compiled nurbforge._core must be importable
    doxygen docs/Doxyfile             # optional: C++ reference (XML in build/docs/xml)
    sphinx-build -b html docs build/sphinx/html

Without Doxygen XML the C++ reference page shows a placeholder note instead of failing.
"""
from __future__ import annotations

import importlib
from pathlib import Path

from docutils import nodes
from docutils.parsers.rst import Directive

DOCS = Path(__file__).resolve().parent
ROOT = DOCS.parent
DOXYGEN_XML = ROOT / "build" / "docs" / "xml"

# -- project information ----------------------------------------------------------------------
project = "Nurbforge"
author = "PhotonPower"
copyright = "PhotonPower, MIT License"

try:
    nurbforge = importlib.import_module("nurbforge")
    release = nurbforge.__version__
    HAVE_PYTHON_PACKAGE = True
except ImportError:  # package not built/installed: autodoc pages will be empty
    import tomllib

    release = tomllib.loads((ROOT / "pyproject.toml").read_text("utf-8"))["project"]["version"]
    HAVE_PYTHON_PACKAGE = False
version = release

# -- general configuration --------------------------------------------------------------------
extensions = [
    "myst_parser",
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.mathjax",
    "sphinx.ext.viewcode",
]

HAVE_DOXYGEN_XML = (DOXYGEN_XML / "index.xml").is_file()
if HAVE_DOXYGEN_XML:
    extensions.append("breathe")
    breathe_projects = {"nurbforge": str(DOXYGEN_XML)}
    breathe_default_project = "nurbforge"
    breathe_default_members = ("members",)

root_doc = "index"
source_suffix = {".md": "markdown", ".rst": "restructuredtext"}
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

myst_enable_extensions = ["colon_fence", "deflist", "dollarmath", "amsmath"]
myst_heading_anchors = 3

autodoc_default_options = {"members": True, "undoc-members": True, "show-inheritance": True}
autodoc_member_order = "bysource"
autodoc_typehints = "description"
# gmsh is an optional runtime dependency (and needs system OpenGL libraries); mock it so the
# nurbforge.mesh docstrings can be imported anywhere.
autodoc_mock_imports = ["gmsh"]

# Pages link to files outside docs/ (README, sources, ADR folder); those links only resolve on
# GitHub. Keep the warnings visible but quiet the noisiest categories.
suppress_warnings = ["myst.header"]

# -- HTML output ------------------------------------------------------------------------------
html_theme = "furo"
html_title = f"Nurbforge {release}"
html_theme_options = {
    "source_repository": "https://github.com/PhotonPower/Nurbforge/",
    "source_branch": "main",
    "source_directory": "docs/",
}


# -- fallbacks when optional inputs are missing ------------------------------------------------
class _Placeholder(Directive):
    """Stands in for Breathe directives when no Doxygen XML is available."""

    has_content = True
    optional_arguments = 1
    final_argument_whitespace = True
    option_spec = None  # accept and ignore any option

    def run(self):
        para = nodes.paragraph(
            text="C++ reference not generated: run `doxygen docs/Doxyfile` from the repository "
            "root before sphinx-build (see CONTRIBUTING.md)."
        )
        return [nodes.note("", para)]


class _AnyOptions(dict):
    def __getitem__(self, key):
        return lambda value: value

    def __contains__(self, key):
        return True


_Placeholder.option_spec = _AnyOptions()


def setup(app):
    if not HAVE_DOXYGEN_XML:
        for name in ("doxygenindex", "doxygennamespace", "doxygenclass", "doxygenstruct",
                     "doxygenfunction", "doxygenfile", "doxygengroup", "doxygenenum"):
            app.add_directive(name, _Placeholder)
