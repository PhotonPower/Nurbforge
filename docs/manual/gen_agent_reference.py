"""Generate the compact API reference for language models and coding agents.

Writes ``python/nurbforge/REFERENCE.md``, which ships inside the wheel, is printed by
``python -m nurbforge`` and is returned by ``nurbforge.reference()``. It packs the whole public
API (signatures plus one-paragraph descriptions) and the conventions an agent needs into a single Markdown file, so that an agent
can use an installed ``nurbforge`` without access to the repository or the manual.

The API part comes from the installed package (the same sources as the manual's Appendix A,
see ``gen_reference.py``); the conventions and examples are written below. ``--check`` fails if
the committed file is out of date or an entry has no description.

    python docs/manual/gen_agent_reference.py [--check]
"""
from __future__ import annotations

import argparse
import importlib
import inspect
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT = ROOT / "python" / "nurbforge" / "REFERENCE.md"
sys.path.insert(0, str(HERE))
import gen_reference as g
import reference_notes as notes

missing: list[str] = []

HEADER = """\
# nurbforge {version} — API reference for agents

Headless NURBS CAD kernel (C++20 core, Python bindings): exact curves and surfaces, trimmed
faces, closed shells (B-rep solids), Booleans, fillets, sweeps, fitting, optical surfaces,
planar slicing/hatching, native STEP read/write, gmsh meshing. `import nurbforge as nf`.
This file is the complete public API. Print it with `python -m nurbforge`;
`nf.reference()` returns it as a string.

## Conventions

- Points/vectors: any 3-sequence (tuple, list, 1-D array) or `nf.Vec3`. Results are `Vec3`
  (`v.x`, `x, y, z = v`, `v.to_tuple()`, arithmetic).
- Angles in **radians** (the LLM tools of `nurbforge.llm` use degrees). No fixed length unit;
  STEP I/O declares/converts units (default mm).
- Every constructed curve/surface has parameter domain [0, 1]; the parameter is not arc length.
  Call objects to evaluate: `curve(u)`, `surface(u, v)`.
- Surface control net is row-major: `index = i * num_v + j`, `i` along u.
- `revolve`: u = circumferential, v = profile. `extrude`/`ruled`: u = profile, v = sweep.
  `loft`: u along the sections, v across them.
- `make_circle(center, normal, r)` runs counter-clockwise about `normal`; its start point is
  chosen by the kernel — use `make_arc`/`make_circle_3pt` to control it.
- Trim loops live in the (u, v) space of the surface: outer counter-clockwise, holes clockwise.
- Solids are closed `Shell`s of `TrimmedSurface` faces with outward normals. Faces carry string
  labels (box: bottom/top/front/back/left/right; cylinder: side/bottom/top; ...). Booleans and
  solid operations keep labels; labels become physical groups when meshing.
- Planar sections: loops with the material on the left (outer CCW, holes CW in plane coords).
- All geometry objects are immutable: methods return new objects. They pickle and deepcopy.
  Heavy functions release the GIL; use threads or processes freely.
- `Tolerance(length=1e-10, parametric=1e-10, angular=1e-9)` is the default `Tolerance()`;
  functions taking `tolerance: float` (Booleans, shells, STEP) default to 1e-6 model units.
- Errors: invalid input or failed geometry raises `nf.GeometryError` (a `ValueError`); a failed
  operation leaves its inputs unchanged.
- Optional extras: `pip install "nurbforge[mesh]"` (gmsh, for `nurbforge.mesh`),
  `"nurbforge[llm]"` (Ollama client, for `nurbforge.llm`).

## Typical workflow

```python
import math
import nurbforge as nf

arc = nf.make_arc((0, 0, 0), (1, 0, 0), (0, 0, 1), 10.0, 0.0, math.pi / 2)
dome = nf.revolve(arc, (0, 0, 0), (0, 0, 1), 2 * math.pi)     # NurbsSurface
print(nf.length(arc), nf.area(dome), dome(0.25, 0.5))

box = nf.make_box((-20, -20, 0), (40, 40, 10))                  # closed Shell
bore = nf.make_cylinder((0, 0, -5), (0, 0, 1), 6.0, 20.0)
part = nf.boolean(box, bore, nf.BooleanOp.Difference)
print(part.is_closed, part.face_labels, part.mass_properties().volume)
print(nf.point_in_shell(part, (10, 10, 5)))

sec = nf.section(part, nf.Plane((0, 0, 5), (0, 0, 1)))          # PlanarSection: loops
lines, contours = nf.hatch(sec, (1, 0, 0), 0.5)                  # hatch segments + contours
left, right = nf.split(part, nf.Plane((0, 0, 0), (1, 0, 0)))    # two closed shells

nf.write_step("part.step", part, unit="mm")
shells = nf.read_step("part.step")                               # list[Shell]
# meshing (needs gmsh): nurbforge.mesh.bem_mesh(part, size=2.0, path="part.msh")
```

Free-form: build curves (`make_*`, `interpolate_curve`, `approximate_curve`), surfaces
(`extrude`, `revolve`, `loft`, `sweep`, `sweep2`, `coons`, `gordon`, `interpolate_surface`),
trim them (`split`, `TrimmedSurface(surface, outer_loop)`), sew faces with
`Shell.sew(faces)`, close planar openings with `cap`, give thickness with `thicken`. Check
results with `shell.is_closed`, `shell.free_edges()`, `shell.volume()`; methods are shown
below as `.name(...)`, properties as `.name: type`.

## API

Signature notation: `Vec3` = any 3-sequence; `array_like` = NumPy array or nested lists;
arguments after `*` are keyword-only; `(1) … (2) …` are overloads.
"""


# ---------------------------------------------------------------------------- text helpers
def md(text: str) -> str:
    """LaTeX snippets of reference_notes -> Markdown."""
    for cmd in ("code", "apiref", "toolref", "texttt"):
        text = re.sub(r"\\" + cmd + r"\{((?:[^{}]|\{[^{}]*\})*)\}", r"`\1`", text)
    for cmd in ("emph", "textbf", "textit", "NF"):
        text = re.sub(r"\\" + cmd + r"\{([^{}]*)\}", r"\1", text)
    text = re.sub(r"\s*\(?\\[cC]ref\{[^}]*\}\)?", "", text)
    text = text.replace(r"\NF{}", "nurbforge").replace(r"\NF", "nurbforge")
    for a, b in ((r"\_", "_"), (r"\%", "%"), (r"\&", "&"), (r"\#", "#"), ("~", " "),
                 (r"\le", "≤"), (r"\ge", "≥"), (r"\times", "×"), (r"\to", "→"), (r"\pi", "π"),
                 (r"\pm", "±"), (r"\cdot", "·"), (r"\ldots", "…"), (r"\dots", "…"),
                 (r"\infty", "∞"), (r"\textdegree{}", "°"), ("---", "—"), ("--", "–"),
                 (r"\{", "{"), (r"\}", "}")):
        text = text.replace(a, b)
    text = re.sub(r"\\textsuperscript\{([^}]*)\}", r"^\1", text)
    text = text.replace("$", "")
    return " ".join(text.split())


def doc_text(text: str) -> str:
    text = re.sub(r":(?:meth|func|class|mod|attr|data|exc):`~?", "`", text)
    return " ".join(text.replace("``", "`").split())


def describe(qual: str, text: str) -> str:
    extra = notes.NOTES.get(qual)
    if text:
        return doc_text(text) + (" " + md(extra) if extra else "")
    if extra:
        return md(extra)
    missing.append(qual)
    return ""


def clean_sig(sig: str) -> str:
    return sig.replace("nf::", "")


def signatures(fn, name: str) -> tuple[list[str], str]:
    if g.is_pybind(fn):
        sigs, text = g.pybind_doc(fn)
    else:
        try:
            sig = str(inspect.signature(fn))
        except (TypeError, ValueError):
            sig = "(...)"
        sig = re.sub(r"\(self(, )?", "(", sig).replace("'", "")
        sigs = [g.clean(fn.__name__ + sig)]
        text = inspect.getdoc(fn) or ""
    sigs = [clean_sig(re.sub(r"^(\(\d+\) )?__init__", lambda m: (m.group(1) or "") + name, s))
            for s in sigs]
    if name[:1].isupper():  # constructor
        sigs = [re.sub(r" -> None$", "", s) for s in sigs]
    return sigs, text


def entry(prefix: str, qual: str, fn, name: str) -> str:
    sigs, text = signatures(fn, name)
    body = describe(qual, text)
    if len(sigs) == 1:
        head = f"- `{prefix}{sigs[0]}`"
    else:
        head = "- " + "  \n  ".join(f"`{prefix}{s}`" for s in sigs)
    return head + (f" — {body}" if body else "") + "\n"


# ---------------------------------------------------------------------------------- classes
def is_enum(cls) -> bool:
    return hasattr(cls, "__members__") and "pybind11" in type(cls).__name__


def enum_block(qual: str, cls) -> str:
    text = md(notes.NOTES.get(qual, "")) or describe(qual, "")
    vals = []
    for k in cls.__members__:
        note = md(notes.NOTES.get(f"{qual}.{k}", ""))
        vals.append(f"`{k}`" + (f" ({note})" if note else ""))
    return f"- **`{qual}`** — {text} Values: " + "; ".join(vals) + "\n"


def class_block(prefix: str, cls) -> str:
    qual = prefix + cls.__name__
    if is_enum(cls):
        return enum_block(qual, cls)
    if issubclass(cls, BaseException):
        text = (cls.__doc__ or "") if cls.__doc__ != Exception.__doc__ else ""
        return (f"### `{qual}` (exception, subclass of {cls.__mro__[1].__name__})\n\n"
                f"{describe(qual, text)}\n\n")
    out = [f"### `{qual}`\n\n{describe(qual, inspect.cleandoc(cls.__doc__ or ''))}\n\n"]
    members = vars(cls)
    init = members.get("__init__")
    if init is not None and (g.is_pybind(init) or inspect.isfunction(init)):
        out.append(entry("", qual + ".__init__", init, cls.__name__))
    for name in sorted(members, key=lambda n: (n.startswith("_"), n)):
        if name in g.SKIP or (name.startswith("_") and name not in g.OPERATORS):
            continue
        v, mq = members[name], f"{qual}.{name}"
        if isinstance(v, property):
            fdoc = (v.fget.__doc__ or "") if v.fget else ""
            m = re.search(r"->\s*(.+)$", fdoc.strip().splitlines()[0]) if fdoc.strip() else None
            typ = clean_sig(g.clean(m.group(1))) if m else ""
            text = " ".join((v.__doc__ or "").split())
            if not text and fdoc and not fdoc.strip().startswith("("):
                text = fdoc
            rw = " (settable)" if v.fset is not None else ""
            body = describe(mq, text)
            out.append(f"- `.{name}" + (f": {typ}" if typ else "") + f"`{rw}"
                       + (f" — {body}" if body else "") + "\n")
        elif isinstance(v, (staticmethod, classmethod)):
            out.append(entry(f"{cls.__name__}.", mq, v.__func__, name).replace(
                "- `", "- static `", 1))
        elif name in g.OPERATORS:
            sigs, _ = signatures(v, name)
            op = g.OPERATORS[name].replace("\\code{", "`").replace("}", "`")
            out.append("- " + ", ".join(f"`{s}`" for s in sigs) + f" — {op}\n")
        elif callable(v):
            out.append(entry(".", mq, v, name))
        elif not name.startswith("_") and not inspect.ismodule(v):
            out.append(f"- `{name} = {v!r}` — {describe(mq, '')}\n")
    return "".join(out) + "\n"


# ----------------------------------------------------------------------------------- driver
def generate() -> str:
    import nurbforge as nf

    out = [HEADER.format(version=nf.__version__)]
    public = [n for n in nf.__all__ if n != "annotations"]
    classes = [getattr(nf, n) for n in public if inspect.isclass(getattr(nf, n))]
    funcs = {n: getattr(nf, n) for n in public
             if callable(getattr(nf, n)) and not inspect.isclass(getattr(nf, n))}

    out.append("\n### Functions of `nurbforge`\n")
    grouped = set()
    for group, names in notes.GROUPS:
        out.append(f"\n**{group}**\n\n")
        for n in names:
            if n not in funcs:
                missing.append(f"group {group}: {n} is not a public function")
                continue
            grouped.add(n)
            out.append(entry("", n, funcs[n], n))
    for n in sorted(set(funcs) - grouped):
        missing.append(f"function {n} is in no group (reference_notes.GROUPS)")

    out.append("\n## Classes of `nurbforge`\n\n")
    enums = [c for c in classes if is_enum(c)]
    for cls in sorted((c for c in classes if c not in enums), key=lambda c: c.__name__):
        out.append(class_block("", cls))
    out.append("### Enumerations\n\n")
    for cls in sorted(enums, key=lambda c: c.__name__):
        out.append(enum_block(cls.__name__, cls))

    for modname in ("nurbforge.mesh", "nurbforge.llm"):
        mod = importlib.import_module(modname)
        out.append(f"\n## `{modname}`\n\n{doc_text(mod.__doc__ or '')}\n\n")
        short = modname.split(".")[-1]
        for n in notes.SUBMODULES[modname]:
            obj = mod
            for part in n.split("."):
                obj = getattr(obj, part)
            qual = f"{short}.{n}"
            if inspect.isclass(obj):
                prefix = qual[: -len(obj.__name__)]
                out.append(class_block(prefix, obj))
            else:
                out.append(entry(f"{modname}.{n.rsplit('.', 1)[0]}." if "." in n else
                                 f"{modname}.", qual, obj, obj.__name__))
        out.append("\n")

    out.append(llm_tools())
    return re.sub(r"\n{3,}", "\n\n", "".join(out)).rstrip() + "\n"


def llm_tools() -> str:
    from nurbforge.llm.tools import TOOLS

    out = ["## LLM tools (`nurbforge.llm.tools.TOOLS`)\n\n",
           ("JSON tools for an LLM agent working on a `Workspace` of named objects; call them "
            "with `nurbforge.llm.tools.call_tool(ws, name, arguments)`. Angles in degrees. "
            "`name?` = optional parameter.\n\n")]
    for t in TOOLS:
        props = t.parameters.get("properties", {})
        req = set(t.parameters.get("required", []))
        args = ", ".join(f"{p}{'' if p in req else '?'}: {g._type(s)}" for p, s in props.items())
        out.append(f"- `{t.name}({args})` — {' '.join(t.description.split())}\n")
    return "".join(out)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--check", action="store_true",
                    help="fail if REFERENCE.md is out of date or entries lack a description")
    ap.add_argument("--out", type=Path, default=OUT)
    args = ap.parse_args()
    text = generate()
    if missing:
        print(f"{len(missing)} entries without description:", file=sys.stderr)
        for m in missing:
            print("  " + m, file=sys.stderr)
        if args.check:
            return 1
    if args.check:
        old = args.out.read_text(encoding="utf-8") if args.out.exists() else ""
        if old != text:
            import difflib

            diff = difflib.unified_diff(old.splitlines(), text.splitlines(), "committed",
                                        "generated", n=0, lineterm="")
            print(f"{args.out} is out of date; run docs/manual/gen_agent_reference.py\n" +
                  "\n".join(list(diff)[:60]), file=sys.stderr)
            return 1
        print(f"{args.out} is up to date")
        return 0
    args.out.write_text(text, encoding="utf-8")
    print(f"wrote {args.out} ({len(text)} characters, ~{len(text) // 4} tokens)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
