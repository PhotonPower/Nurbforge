# Nurbforge documentation

| Document | Content |
|---|---|
| [architecture.md](architecture.md) | Layers, modules, dependencies, error handling |
| [math/nurbs.md](math/nurbs.md) | NURBS definitions, conventions, implemented algorithms |
| [api/python.md](api/python.md) | Python API overview |
| [api/cpp.md](api/cpp.md) | C++ API overview, Doxygen |
| [api/python-reference.md](api/python-reference.md) | Python API reference (autodoc, rendered on the docs site) |
| [api/cpp-reference.md](api/cpp-reference.md) | C++ API reference (Doxygen + Breathe, rendered on the docs site) |
| [optics.md](optics.md) | Optical surfaces: conics, aspheres, Zernike, XY polynomials |
| [meshing.md](meshing.md) | gmsh integration |
| [llm.md](llm.md) | Ollama interface and safety model |
| [build/almalinux10.md](build/almalinux10.md) | Build on AlmaLinux 10 |
| [build/windows11.md](build/windows11.md) | Build on Windows 11 |
| [roadmap.md](roadmap.md) | Planned features and milestones |
| [adr/README.md](adr/README.md) | Architecture decision records |

Standalone Doxygen HTML: run `doxygen docs/Doxyfile` from the repository root; output goes to
`build/docs/html` (XML for the Sphinx site in `build/docs/xml`).

The rendered documentation site (Sphinx, built from these Markdown files plus the generated API
references) is published at <https://photonpower.github.io/Nurbforge/>. To build it locally see
[CONTRIBUTING.md](https://github.com/PhotonPower/Nurbforge/blob/main/CONTRIBUTING.md#documentation).

<!-- Navigation tree for the Sphinx site (docs/conf.py); GitHub shows it as a plain code block. -->

```{toctree}
:hidden:
:caption: Guide

architecture
math/nurbs
optics
meshing
llm
build/almalinux10
build/windows11
```

```{toctree}
:hidden:
:caption: API

api/python
api/python-reference
api/cpp
api/cpp-reference
```

```{toctree}
:hidden:
:caption: Project
:glob:

roadmap
adr/README
adr/0*
```
