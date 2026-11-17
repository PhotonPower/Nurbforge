#pragma once
/// @file solid_ops.hpp
/// Operations that build or modify solids for simulation: capping openings, thickening, ...

#include <functional>
#include <string>

#include "nurbforge/modeling/sweep.hpp"
#include "nurbforge/topology/shell.hpp"

namespace nf {

/// Closes the planar openings of a shell: the free edges are chained into closed loops, each loop
/// must lie in a plane (within `tolerance`), and loops in the same plane become one planar face
/// (nested loops: holes, e.g. the annular end of a tube; islands inside holes are new faces).
/// The faces are planar rectangles trimmed exactly by the opening's edges (the trim curves are
/// the edge curves mapped affinely into the rectangle's parameter space). They get the label
/// `label`; the result is sewn and oriented outwards. A closed shell is returned unchanged.
/// Throws GeometryError if a loop does not close or is not planar.
[[nodiscard]] Shell cap(const Shell& shell, double tolerance = 1e-6,
                        const std::string& label = "cap");

/// Gives an open shell (or a single face, sewn into a shell) a wall thickness: every face is
/// offset along its normal S_u x S_v by `thickness` (negative: against the normal; with
/// `symmetric` by +-thickness/2 on both sides), keeping its trim loops, and every free edge is
/// closed by a ruled side wall between the edge and its offset. The result is sewn and oriented
/// outwards; the offset faces keep their labels with the suffix "_offset" (unlabelled faces:
/// "offset"), the side walls get
/// `wall_label`. The faces must join smoothly (tangent-continuous) along their shared edges, so
/// that their offsets meet; at a crease the result is not closed and GeometryError is raised, as
/// it is where |thickness| reaches a radius of curvature (offset_surface).
[[nodiscard]] Shell thicken(const Shell& shell, double thickness, bool symmetric = false,
                            double tolerance = 1e-6, const std::string& wall_label = "wall");

/// Solid pipe along `path`: the pipe surface (pipe_surface, radius(t) for t in [0, 1]) closed by
/// planar discs at both ends (cap). Faces are labelled "wall", "inlet" (path start) and "outlet"
/// (path end); a closed path gives a ring without end faces. `tolerance` is the deviation of
/// the wall from the ideal pipe and the sewing tolerance is 100 times larger.
[[nodiscard]] Shell make_pipe(const NurbsCurve& path, const std::function<double(double)>& radius,
                              double tolerance = 1e-6);
/// Pipe with a radius varying linearly from `start_radius` to `end_radius` (negative end radius:
/// constant).
[[nodiscard]] Shell make_pipe(const NurbsCurve& path, double start_radius, double end_radius = -1.0,
                              double tolerance = 1e-6);

}  // namespace nf
