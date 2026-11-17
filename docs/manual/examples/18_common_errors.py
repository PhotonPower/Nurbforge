"""Typical mistakes and the messages they produce."""
import math

import nurbforge as nf


def attempt(label, fn):
    try:
        fn()
        print(f"{label}: ok")
    except nf.GeometryError as exc:
        print(f"{label}: GeometryError: {exc}")
    except (TypeError, ValueError) as exc:
        print(f"{label}: {type(exc).__name__}: {str(exc).splitlines()[0]}")


# wrong number of knots for the control points
attempt("knot vector", lambda: nf.NurbsCurve(3, [0, 0, 1, 1], [(0, 0, 0), (1, 0, 0)]))
# a zero-length axis
attempt("zero axis", lambda: nf.make_cylinder((0, 0, 0), (0, 0, 0), 1.0, 2.0))
# an offset larger than the radius of curvature on the concave side
circle = nf.make_circle((0, 0, 0), (0, 0, 1), 1.0)
attempt("offset > radius", lambda: nf.offset_curve(circle, 1.5))
# a point given with two coordinates
attempt("2D point", lambda: nf.make_line((0, 0), (1, 1)))
# angles are radians: 90 means 90 rad, i.e. more than 14 turns
attempt("angle in degrees", lambda: nf.revolve(nf.make_line((1, 0, 0), (1, 0, 1)),
                                               (0, 0, 0), (0, 0, 1), 90))
attempt("angle in radians", lambda: nf.revolve(nf.make_line((1, 0, 0), (1, 0, 1)),
                                               (0, 0, 0), (0, 0, 1), math.radians(90)))
# Booleans need closed shells
box = nf.make_box((0, 0, 0), (1, 1, 1))
open_shell = nf.Shell.sew(list(box.faces)[:-1])
print("open shell closed:", open_shell.is_closed, " free edges:", len(open_shell.free_edges()))
attempt("boolean with open shell",
        lambda: nf.boolean(open_shell, nf.make_sphere((1, 1, 1), 0.5), nf.BooleanOp.Union))
