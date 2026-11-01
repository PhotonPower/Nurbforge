"""Interactive command line: ``nurbforge-llm [--model M] [--host H] [--out file.msh]``."""
from __future__ import annotations

import argparse

from .agent import GeometryAgent
from .ollama_client import DEFAULT_HOST, DEFAULT_MODEL, OllamaClient


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Build Nurbforge geometry by talking to a local LLM.")
    ap.add_argument("--model", default=DEFAULT_MODEL)
    ap.add_argument("--host", default=DEFAULT_HOST)
    ap.add_argument("--out", help="mesh all surfaces with gmsh and write to this file on exit")
    ap.add_argument("--mesh-size", type=float, default=None)
    ap.add_argument("--verify-rounds", type=int, default=1,
                    help="verification rounds after each answer (0 disables)")
    args = ap.parse_args(argv)

    agent = GeometryAgent(OllamaClient(args.host, args.model), verify_rounds=args.verify_rounds)
    print(f"Nurbforge LLM shell ({args.model} @ {args.host}). Empty line to quit.")
    while True:
        try:
            prompt = input("> ").strip()
        except EOFError:
            break
        if not prompt:
            break
        print(agent.ask(prompt))
        for name, a, res in agent.log[-10:]:
            print(f"  [{name}] {a} -> {res}")
        agent.log.clear()

    if args.out:
        from .. import NurbsSurface
        from ..mesh import mesh_surfaces
        surfaces = [o for o in agent.workspace.objects.values() if isinstance(o, NurbsSurface)]
        mesh_surfaces(surfaces, args.out, size=args.mesh_size)
        print(f"wrote {len(surfaces)} surface(s) to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
