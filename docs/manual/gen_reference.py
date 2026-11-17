"""Generate the command reference of the user manual from the installed ``nurbforge`` package.

Writes ``generated/reference.tex`` (Appendix A: every public class, method, property, enum and
function of ``nurbforge``, ``nurbforge.mesh`` and ``nurbforge.llm``), ``generated/llm-tools.tex``
(Appendix B: every LLM tool with its parameters) and ``generated/version.tex``.

Descriptions come from the docstrings; entries without one take their text from
``reference_notes.py``. ``--check`` fails if any public entry has no description or a function is
not assigned to a group, so the reference stays complete when the API grows.

    python docs/manual/gen_reference.py [--check] [--out DIR]
"""
from __future__ import annotations

import argparse
import importlib
import inspect
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import reference_notes as notes

OPERATORS = {
    "__call__": "Call operator.", "__add__": "Operator \\code{a + b}.",
    "__sub__": "Operator \\code{a - b}.", "__mul__": "Operator \\code{a * b}.",
    "__rmul__": "Operator \\code{b * a}.", "__truediv__": "Operator \\code{a / b}.",
    "__neg__": "Operator \\code{-a}.", "__getitem__": "Indexing \\code{a[i]}.",
    "__len__": "\\code{len(a)}.", "__iter__": "Iteration \\code{for x in a}.",
    "__eq__": "Comparison \\code{a == b}.", "__ne__": "Comparison \\code{a != b}.",
}
SKIP = {"__init__", "__repr__", "__str__", "__hash__", "__int__", "__index__", "__setstate__",
        "__getstate__", "__reduce__", "__new__", "__module__", "__doc__", "__dict__",
        "__weakref__", "__members__", "__entries", "_pybind11_conduit_v1_", "__annotations__",
        "__dataclass_fields__", "__dataclass_params__", "__match_args__", "__firstlineno__",
        "__static_attributes__", "__enter__", "__exit__", "__qualname__", "__orig_bases__",
        "__parameters__", "__abstractmethods__", "__slots__", "__class_getitem__",
        "__annotate_func__", "__annotations_cache__"}

missing: list[str] = []


# ------------------------------------------------------------------------------------ text
def tex(text: str) -> str:
    """Docstring text -> LaTeX: escape specials, `code` and ``code`` -> \\code{}."""
    text = re.sub(r":(?:meth|func|class|mod|attr|data|exc):`~?", "`", text)
    out, pos = [], 0
    for m in re.finditer(r"``(.+?)``|`([^`]+)`", text):
        out.append(_esc(text[pos:m.start()]))
        out.append(_code(m.group(1) or m.group(2)))
        pos = m.end()
    out.append(_esc(text[pos:]))
    return "".join(out)


def _esc(s: str) -> str:
    rep = {"\\": r"\textbackslash{}", "{": r"\{", "}": r"\}", "$": r"\$", "&": r"\&",
           "#": r"\#", "%": r"\%", "_": r"\_", "^": r"\^{}", "~": r"\~{}",
           "<": r"\textless{}", ">": r"\textgreater{}"}
    s = "".join(rep.get(c, c) for c in s)
    for a, b in {"→": r"$\to$", "×": r"$\times$", "²": r"\textsuperscript{2}",
                 "³": r"\textsuperscript{3}", "π": r"$\pi$", "≤": r"$\le$", "≥": r"$\ge$",
                 "·": r"$\cdot$", "µ": r"\textmu{}", "–": "--", "—": "---", "°": r"\textdegree{}",
                 "≈": r"$\approx$", "±": r"$\pm$", "√": r"$\surd$", "∈": r"$\in$",
                 "⁻": "-", "¹": r"\textsuperscript{1}", "ρ": r"$\rho$", "θ": r"$\theta$",
                 "α": r"$\alpha$", "Δ": r"$\Delta$", "∞": r"$\infty$", "≠": r"$\ne$"}.items():
        s = s.replace(a, b)
    return s


def _code(s: str) -> str:
    # \texttt with escaped characters (\code{} would choke on unbalanced braces)
    return r"\texttt{" + _esc(s) + "}"


def lst(code: str) -> str:
    code = code.replace("\\end{lstlisting}", "\\end {lstlisting}")
    return "\\begin{lstlisting}[style=sig]\n" + code + "\n\\end{lstlisting}\n"


# ------------------------------------------------------------------------------ signatures
_TYPES = [
    (r"typing\.SupportsFloat \| typing\.SupportsIndex", "float"),
    (r"typing\.SupportsInt \| typing\.SupportsIndex", "int"),
    (r'typing\.Annotated\[numpy\.typing\.ArrayLike, numpy\.float64\]', "array_like"),
    (r'typing\.Annotated\[numpy\.typing\.NDArray\[numpy\.float64\], [^\]]*\]', "ndarray"),
    (r"numpy\.typing\.NDArray\[numpy\.float64\]", "ndarray"),
    (r'typing\.Annotated\[(list|Sequence)\[([^\]]*)\], "FixedSize\((\d+)\)"\]', r"\1[\2] (\3)"),
    (r'typing\.Annotated\[(.*?), "FixedSize\((\d+)\)"\]', r"\1 (\2)"),
    (r"nurbforge\._core\.", ""), (r"collections\.abc\.", ""), (r"typing\.", ""),
    (r"numpy\.ndarray", "ndarray"), (r"pybind11_builtins\.", ""),
    (r"Tolerance\(length=0\.000000, parametric=0\.000000, angular=0\.000000\)", "Tolerance()"),
    (r"<(\w+)\.(\w+): \d+>", r"\1.\2"),
    (r"<(\w+) object at 0x[0-9a-f]+>", r"\1()"),
    (r"Vec3\((-?[\d.]+), (-?[\d.]+), (-?[\d.]+)\)",
     lambda m: "(" + ", ".join(_num(x) for x in m.groups()) + ")"),
]


def _num(x: str) -> str:
    f = float(x)
    return str(int(f)) if f == int(f) else repr(f)


def clean(sig: str) -> str:
    for pat, rep in _TYPES:
        sig = re.sub(pat, rep, sig)
    return re.sub(r"\(self(: [\w.]+)?(, )?", "(", sig, count=1)


_SIG = re.compile(r"^\s*(?:\d+\.\s+)?([A-Za-z_]\w*)\((.*)\)(?:\s*->\s*(.+))?\s*$")


def pybind_doc(obj) -> tuple[list[str], str]:
    """Signatures and description of a pybind11 function from its docstring. Overloads keep
    their own descriptions, numbered like the signatures."""
    doc = inspect.cleandoc(getattr(obj, "__doc__", None) or "")
    name = getattr(obj, "__name__", "")
    blocks: list[tuple[str, list[str]]] = []
    lead: list[str] = []
    for line in doc.splitlines():
        if line.strip() == "Overloaded function.":
            continue
        m = _SIG.match(line)
        if m and m.group(1) == name:
            if m.group(2) == "*args, **kwargs":
                continue
            blocks.append((clean(re.sub(r"^\s*\d+\.\s+", "", line).strip()), []))
        elif blocks:
            blocks[-1][1].append(line)
        else:
            lead.append(line)
    sigs = [b[0] for b in blocks]
    texts = [" ".join(" ".join(b[1]).split()) for b in blocks]
    if len(blocks) > 1 and any(texts):
        text = " ".join(f"({i}) {t}" for i, t in enumerate(texts, 1) if t)
        sigs = [f"({i}) {s}" for i, s in enumerate(sigs, 1)]
    else:
        text = " ".join(t for t in texts if t)
    text = " ".join((" ".join(lead) + " " + text).split())
    return sigs, text


def is_pybind(obj) -> bool:
    return type(obj).__name__ in ("builtin_function_or_method", "instancemethod") or \
        "pybind11" in type(obj).__name__


def describe(qual: str, text: str) -> str:
    if text:
        extra = notes.NOTES.get(qual)
        return tex(text) + (" " + extra if extra else "")
    if qual in notes.NOTES:
        return notes.NOTES[qual]
    missing.append(qual)
    return ""


# --------------------------------------------------------------------------------- entries
def function_entry(qual: str, fn, kind: str = "function") -> str:
    if is_pybind(fn):
        sigs, text = pybind_doc(fn)
        name = fn.__name__
        sigs = [s if not s.startswith("__init__") else s.replace("__init__", qual.split(".")[-2]
                if qual.endswith("__init__") else name, 1) for s in sigs]
    else:
        try:
            sig = str(inspect.signature(fn))
        except (TypeError, ValueError):
            sig = "(...)"
        sig = re.sub(r"\(self(, )?", "(", sig)
        sigs = [clean(fn.__name__ + sig)]
        text = " ".join((inspect.getdoc(fn) or "").split())
    shown = qual if not qual.endswith(".__init__") else qual[: -len(".__init__")] + "()"
    if qual.endswith(".__init__"):
        cls = shown.split(".")[-1]
        sigs = [re.sub(r"^__init__", cls, s) for s in sigs]
        sigs = [re.sub(r"\) -> None$", ")", s) for s in sigs]
    body = describe(qual, text)
    return (f"\\begin{{refentry}}{{{shown}}}{{{kind}}}\n" + lst("\n".join(sigs)) + body +
            "\n\\end{refentry}\n")


def property_entry(qual: str, prop) -> str:
    fdoc = (prop.fget.__doc__ or "") if prop.fget else ""
    m = re.search(r"->\s*(.+)$", fdoc.strip().splitlines()[0]) if fdoc.strip() else None
    typ = clean(m.group(1)) if m else ""
    text = " ".join((prop.__doc__ or "").split())
    if not text and fdoc and not fdoc.strip().startswith("("):
        text = " ".join(fdoc.split())
    ro = "" if prop.fset is not None else " (read-only)"
    sig = qual.split(".")[-1] + (f": {typ}" if typ else "")
    return (f"\\begin{{refentry}}{{{qual}}}{{property{ro}}}\n" + lst(sig) +
            describe(qual, text) + "\n\\end{refentry}\n")


def enum_entry(qual: str, cls) -> str:
    vals = ", ".join(cls.__members__)
    text = notes.NOTES.get(qual) or describe(qual, "")
    items = "".join(f"\\item[\\texttt{{{_esc(k)}}}] {notes.NOTES.get(qual + '.' + k, '')}\n"
                    for k in cls.__members__)
    return (f"\\begin{{refentry}}{{{qual}}}{{enum}}\n" + lst(f"{qual.split('.')[-1]}.{{{vals}}}")
            + text + "\n\\begin{description}[style=nextline,leftmargin=2em,itemsep=0pt]\n" +
            items + "\\end{description}\n\\end{refentry}\n")


def class_section(prefix: str, cls) -> str:
    qual = prefix + cls.__name__
    if hasattr(cls, "__members__") and "pybind11" in type(cls).__name__:
        return enum_entry(qual, cls)
    parts = [(f"\\subsection{{\\texorpdfstring{{\\texttt{{{_esc(cls.__name__)}}}}}"
              f"{{{cls.__name__}}}}}\n")]
    if issubclass(cls, BaseException):
        text = " ".join((cls.__doc__ or "").split()) if cls.__doc__ and \
            cls.__doc__ != Exception.__doc__ else ""
        base = cls.__mro__[1].__name__
        parts.append(f"\\begin{{refentry}}{{{qual}}}{{exception}}\n"
                     + lst(f"class {cls.__name__}({base})") + describe(qual, text) +
                     "\n\\end{refentry}\n")
        return "".join(parts)
    doc = inspect.cleandoc(cls.__doc__ or "")
    if is_pybind(getattr(cls, "__init__", None)) or "pybind11" in type(cls).__name__:
        text = " ".join(doc.split())
    else:
        text = " ".join(doc.split())
    parts.append(f"\\begin{{refentry}}{{{qual}}}{{class}}\n" + lst(f"class {cls.__name__}") +
                 describe(qual, text) + "\n\\end{refentry}\n")
    members = vars(cls)
    init = members.get("__init__")
    if init is not None and (is_pybind(init) or inspect.isfunction(init)):
        parts.append(function_entry(qual + ".__init__", init, "constructor"))
    for name in sorted(members, key=lambda n: (n.startswith("_"), n)):
        if name in SKIP or (name.startswith("_") and name not in OPERATORS):
            continue
        v = members[name]
        mq = f"{qual}.{name}"
        if isinstance(v, property):
            parts.append(property_entry(mq, v))
        elif isinstance(v, (staticmethod, classmethod)):
            parts.append(function_entry(mq, v.__func__, "static method"))
        elif name in OPERATORS:
            fn = v
            sigs, _ = pybind_doc(fn) if is_pybind(fn) else ([clean(name + str(
                inspect.signature(fn)))], "")
            parts.append(f"\\begin{{refentry}}{{{mq}}}{{operator}}\n" + lst("\n".join(sigs))
                         + OPERATORS[name] + "\n\\end{refentry}\n")
        elif callable(v):
            parts.append(function_entry(mq, v, "method"))
        elif not name.startswith("_") and not inspect.ismodule(v):
            parts.append(f"\\begin{{refentry}}{{{mq}}}{{attribute}}\n" + lst(
                f"{name} = {v!r}") + describe(mq, "") + "\n\\end{refentry}\n")
    return "".join(parts)


# ---------------------------------------------------------------------------------- driver
def reference() -> str:
    import nurbforge as nf

    out = [r"""\chapter{Command reference}\label{app:reference}
This appendix lists every public class, method, property, enumeration and function of the
Python package: \code{nurbforge} (imported as \code{nf}), \code{nurbforge.mesh} and
\code{nurbforge.llm}. It is generated from the installed package (\code{gen_reference.py}), so
the signatures are exactly those of version \NFVersion. Points and vectors written as
\code{Vec3} also accept tuples and lists of three numbers; \code{array_like} accepts NumPy
arrays and nested lists; \code{Tolerance()} denotes the default tolerances
(\cref{app:tolerances}). Arguments after \code{*} are keyword-only.
"""]
    public = [n for n in nf.__all__ if n not in ("annotations",)] if hasattr(nf, "__all__") \
        else [n for n in dir(nf) if not n.startswith("_")]
    classes = [getattr(nf, n) for n in public if inspect.isclass(getattr(nf, n))]
    funcs = {n: getattr(nf, n) for n in public
             if callable(getattr(nf, n)) and not inspect.isclass(getattr(nf, n))}

    out.append("\\section{Functions of \\texttt{nurbforge}}\\label{app:ref:functions}\n")
    grouped = set()
    for group, names in notes.GROUPS:
        out.append(f"\\subsection{{{group}}}\n")
        for n in names:
            if n not in funcs:
                missing.append(f"group {group}: {n} is not a public function")
                continue
            grouped.add(n)
            out.append(function_entry(n, funcs[n]))
    for n in sorted(set(funcs) - grouped):
        missing.append(f"function {n} is in no group (reference_notes.GROUPS)")

    out.append("\\section{Classes and enumerations of \\texttt{nurbforge}}"
               "\\label{app:ref:classes}\n")
    enums = [c for c in classes if hasattr(c, "__members__") and "pybind11" in type(c).__name__]
    for cls in sorted((c for c in classes if c not in enums), key=lambda c: c.__name__):
        out.append(class_section("", cls))
    out.append("\\subsection{Enumerations}\n")
    for cls in sorted(enums, key=lambda c: c.__name__):
        out.append(enum_entry(cls.__name__, cls))

    for modname, title, label in (("nurbforge.mesh", "Meshing: \\texttt{nurbforge.mesh}", "mesh"),
                                  ("nurbforge.llm", "LLM interface: \\texttt{nurbforge.llm}",
                                   "llm")):
        mod = importlib.import_module(modname)
        out.append(f"\\section{{{title}}}\\label{{app:ref:{label}}}\n")
        out.append(tex(" ".join((mod.__doc__ or "").split())) + "\n")
        names = notes.SUBMODULES[modname]
        for n in names:
            obj = mod
            for part in n.split("."):
                obj = getattr(obj, part)
            if inspect.isclass(obj):
                out.append(class_section(label + "." + ".".join(n.split(".")[:-1]) +
                                         ("." if "." in n else ""), obj)
                           .replace(f"{{{label}..", f"{{{label}."))
            else:
                out.append(function_entry(f"{label}.{n}", obj))
        exported = {k for k, v in vars(mod).items() if not k.startswith("_") and
                    (inspect.isclass(v) or inspect.isfunction(v))
                    and getattr(v, "__module__", "").startswith(modname)}
        for n in sorted(exported - {x.split(".")[-1] for x in names} - {"TOOLS"}):
            missing.append(f"{modname}.{n} is not listed in reference_notes.SUBMODULES")
    return "".join(out)


def llm_tools() -> str:
    from nurbforge.llm.tools import TOOLS

    out = [r"""\chapter{LLM tool reference}\label{app:llmtools}
The tools a language model may call through \code{nurbforge.llm} (\cref{ch:llm}), in the order
in which they are offered to the model. Every tool works on named objects of the
\code{Workspace}; \emph{name} arguments refer to objects stored there. Angles are in degrees.
Points and vectors are arrays of three numbers. The descriptions are the ones the model sees.
"""]
    for t in TOOLS:
        props = t.parameters.get("properties", {})
        req = set(t.parameters.get("required", []))
        out.append(f"\\begin{{toolentry}}{{{t.name}}}\n")
        out.append(tex(t.description) + "\n")
        if props:
            out.append("\n\\begin{tabularx}{\\linewidth}{@{}>{\\ttfamily}l l l X@{}}\n"
                       "\\toprule\\textrm{Parameter} & Type & Required & Description\\\\\\midrule\n")
            for p, spec in props.items():
                out.append(f"{_esc(p)} & {_esc(_type(spec))} & {'yes' if p in req else 'no'} & "
                           f"{tex(spec.get('description', ''))}\\\\\n")
            out.append("\\bottomrule\n\\end{tabularx}\n")
        out.append("\\end{toolentry}\n")
    return "".join(out)


def _type(spec: dict) -> str:
    t = spec.get("type", "any")
    if "enum" in spec:
        return " | ".join(map(str, spec["enum"]))
    if t == "array":
        item = _type(spec.get("items", {}))
        n = spec.get("minItems")
        return f"{item}[{n}]" if n and n == spec.get("maxItems") else f"{item}[]"
    return t if isinstance(t, str) else " | ".join(t)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--check", action="store_true", help="fail on undocumented entries")
    ap.add_argument("--out", type=Path, default=HERE / "generated")
    args = ap.parse_args()
    import nurbforge as nf

    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "reference.tex").write_text(reference(), encoding="utf-8")
    (args.out / "llm-tools.tex").write_text(llm_tools(), encoding="utf-8")
    (args.out / "version.tex").write_text(f"\\newcommand{{\\NFVersion}}{{{nf.__version__}}}\n",
                                          encoding="utf-8")
    if missing:
        print(f"{len(missing)} reference entries without description:", file=sys.stderr)
        for m in missing:
            print("  " + m, file=sys.stderr)
        if args.check:
            return 1
    print(f"wrote {args.out}/reference.tex, llm-tools.tex, version.tex (nurbforge {nf.__version__})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
