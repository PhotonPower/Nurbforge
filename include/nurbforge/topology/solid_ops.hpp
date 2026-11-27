#pragma once
/// @file solid_ops.hpp
/// Operations that build or modify solids for simulation: capping openings, thickening, ...

#include <functional>
#include <string>
#include <vector>

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

/// Defeaturing: removes holes smaller than `max_size` from a closed shell. Every hole loop of a
/// face whose 3D extent (bounding-box diagonal) is below `max_size` is dropped – the face is
/// filled there – and the faces that formed the feature (the bore wall of a through hole; wall
/// and bottom of a blind hole) are removed: after dropping the loops they are attached to the
/// rest only along free edges, and they are deleted repeatedly until the shell is closed again.
/// Throws GeometryError if the result is not closed (e.g. a feature that is not a hole).
[[nodiscard]] Shell remove_holes(const Shell& shell, double max_size, double tolerance = 1e-6);

/// Hollows a smooth closed shell: every face is offset by `thickness` against its outward
/// normal, and the offset shell becomes an inner void (Shell with two components). With
/// `open_labels` the faces carrying these labels are removed and the rest is thickened inwards
/// instead (thicken), giving a cup-like part open at those faces. Like thicken, the faces must
/// join smoothly, so that their offsets meet (GeometryError otherwise).
[[nodiscard]] Shell hollow(const Shell& shell, double thickness,
                           const std::vector<std::string>& open_labels = {},
                           double tolerance = 1e-6);

/// Defeaturing: merges adjacent faces that lie on the same surface and carry the same label
/// into one face – e.g. the coplanar pieces of the top faces after the union of two boxes.
/// Faces qualify if they share an edge and either have the same surface (pieces of one split
/// face) or are planar parallelograms in one plane with the same normal; the merged face is the
/// common surface or a rectangle in that plane, trimmed by the union of the regions (shared
/// boundary pieces cancel, T-junctions are resolved). Trim curves stay exact (affine maps
/// between planar parametrisations). Groups whose boundary does not form valid loops are left
/// as they are. The result is sewn and oriented like the input.
[[nodiscard]] Shell merge_faces(const Shell& shell, double tolerance = 1e-6);

/// Slab: a planar curve given a width and a height. The curve is offset in its plane (normal
/// to `direction`) by `width` (> 0: to the left of the curve seen from the tip of `direction`,
/// < 0: to the right), corners mitred as in ribbon(), and the strip between the curve and its
/// offset is extruded along `direction` (its length is the height) into a closed solid: a
/// wall along an open curve, a ring along a closed one. Faces are labelled "side" (the walls
/// along the curve, its offset and the ends) and "cap" (bottom and top). Throws GeometryError
/// if the curve is not planar normal to the direction or the offset fails (offset_curve).
[[nodiscard]] Shell make_slab(const NurbsCurve& curve, double width, const Vec3& direction,
                              double tolerance = 1e-6);

/// Wire cut: cuts a closed shell along a curve extruded in `direction`, like a wire through
/// foam. The curve is projected along the direction onto a plane; a closed curve cuts out the
/// prism over its interior, an open curve must run right across the solid and divides it
/// into the part to its left (seen from the tip of `direction`) and the part to its right.
/// Returns the parts that are not empty, the left (or inner) one first; the cut faces get the
/// label `cut_label`, all other faces keep theirs.
///
/// Implemented with Booleans against a prism: over the closed curve, or over the open curve
/// closed far away on its left by three straight lines from its end points. For an open curve
/// this loop must not cross itself, so the curve should not curl back across the lines from
/// its ends (a line, an arc, a gentle spline are fine).
[[nodiscard]] std::vector<Shell> wire_cut(const Shell& shell, const NurbsCurve& curve,
                                          const Vec3& direction, double tolerance = 1e-6,
                                          const std::string& cut_label = "cut");

}  // namespace nf
