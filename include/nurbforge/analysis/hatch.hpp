#pragma once

/// @file hatch.hpp
/// Hatching of planar section loops: parallel hatch lines clipped to the region bounded by the
/// loops, and inward contour lines (offsets of the loops). Used to fill the layers of a sliced
/// body with writing paths, e.g. for two-photon lithography (2PP), where a layer is written as
/// contour lines plus parallel hatch lines. All functions are const queries and thread-safe.

#include <cstddef>
#include <cstdint>
#include <vector>

#include "nurbforge/analysis/shell_section.hpp"
#include "nurbforge/core/vec.hpp"
#include "nurbforge/nurbs/curve.hpp"

namespace nf {

struct HatchOptions {
  /// Offset of the line family along the line normal n (see hatch): the lines are
  /// n . p = phase + k * spacing for integer k, measured from the origin of the loop coordinates,
  /// so the same phase gives the same lines in every layer.
  double phase{0.0};
  /// Inward offset distances (> 0) of the contour lines; one closed contour per loop and distance.
  std::vector<double> contours;
  /// Hatch segments shorter than this are dropped.
  double min_length{0.0};
  /// Length tolerance (model units) for contacts: where a loop comes within this distance of a
  /// hatch line without crossing it (tangential contact, an edge lying on the line) the line is
  /// not split. 0: 1e-9 of the model size.
  double tolerance{0.0};
  /// Fit tolerance of the contour curves (see offset_curve). 0: 1e-7 of the model size.
  double contour_tolerance{0.0};
};

/// One hatch segment, inside the region bounded by the loops.
struct HatchSegment {
  Vec3 start, end;            ///< end points (x, y, 0), start before end along the direction
  std::int64_t line{0};       ///< index k of the hatch line n . p = phase + k * spacing
  std::size_t region{0};      ///< index of the outer loop (area > 0) of the region containing it
  std::size_t start_loop{0};  ///< loop on which the segment starts
  std::size_t end_loop{0};    ///< loop on which the segment ends
  [[nodiscard]] double length() const { return distance(start, end); }
};

/// Closed contour line: boundary of the region shrunk by `distance`, oriented like the loops
/// (region on the left). Usually the offset of one loop; where the region splits (a neck narrower
/// than 2 distance) a loop gives several contours, and where a hole comes closer than 2 distance
/// to another loop their contours merge into one.
struct HatchContour {
  std::vector<NurbsCurve> curves;  ///< consecutive pieces (offsets, arcs at reflex corners)
  std::size_t loop{0};             ///< loops.front()
  std::vector<std::size_t> loops;  ///< indices of the loops it was offset from (sorted)
  std::size_t index{0};            ///< index of the distance in HatchOptions::contours
  double distance{0.0};
};

struct Hatching {
  /// Sorted by line index, then along the direction.
  std::vector<HatchSegment> segments;
  /// Sorted by (first) source loop, then by the order of the distances in HatchOptions::contours.
  std::vector<HatchContour> contours;
};

/// Hatching of the region bounded by closed planar loops with lines along `direction`
/// (projected into the xy plane and normalised; z of the loops is ignored) at distance
/// `spacing`.
///
/// Each loop is a chain of curves in the xy plane (end of one = start of the next, last ends at
/// the first), oriented with the region on the left: outer loops counter-clockwise, holes
/// clockwise (as returned by section). Islands inside holes are further counter-clockwise loops.
/// The region is decided by this orientation (winding number), not by even-odd counting.
///
/// The lines are n . p = phase + k spacing with n = (-d.y, d.x). Every line is intersected
/// exactly with the loop curves: the height n . c(t) of each Bezier piece is a rational
/// function with the same weights, so its sign changes are found by subdivision of the
/// control values (single sign change: one transversal root, refined to machine precision).
/// A loop that only touches a line within the tolerance, or runs along it and returns to the
/// same side (an edge on the line), does not split it; where the loop runs along a line and
/// then crosses it, the segment ends where the loop leaves the line to the side n . p > c.
/// Consequently tangent lines give no zero-length segments, and edges lying exactly on a hatch
/// line are not hatched themselves (write them as contours).
///
/// Contours (options.contours) are the boundaries of the region shrunk by each distance d, i.e.
/// of the points of the region at distance >= d from all loops. The raw offset of every loop
/// (offset_curve of its pieces to the left, i.e. into the region, also for holes; circular arcs
/// about reflex corners; parts where d reaches the radius of curvature of a concave piece left
/// out) is split at all its self- and mutual intersections; the fragments at distance >= d from
/// the loops (within the fit tolerance) and inside the region are chained into closed contours.
/// So contours are trimmed at convex corners, split where the region is narrower than 2 d,
/// merged where a hole comes closer than 2 d to another loop, and vanish where the region is
/// thinner than 2 d everywhere (no contour for that loop). Consecutive contour curves meet within
/// max(100 contour_tolerance, 10 x the largest gap between consecutive loop curves).
///
/// Throws GeometryError if spacing <= 0, the direction is zero in the xy plane, or the number of
/// lines exceeds 10^7.
[[nodiscard]] Hatching hatch(const std::vector<std::vector<NurbsCurve>>& loops,
                             const Vec3& direction, double spacing,
                             const HatchOptions& options = {});

/// Hatching of section loops (their curves in the plane coordinates (s, t, 0)); same as above.
[[nodiscard]] Hatching hatch(const std::vector<SectionLoop>& loops, const Vec3& direction,
                             double spacing, const HatchOptions& options = {});

}  // namespace nf
