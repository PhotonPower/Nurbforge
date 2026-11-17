#pragma once
/// @file region.hpp
/// Planar regions bounded by closed curves, and Boolean operations on them (2D profiles,
/// cross-sections, domains for 2D meshing).

#include <vector>

#include "nurbforge/nurbs/curve.hpp"
#include "nurbforge/topology/boolean.hpp"

namespace nf {

/// A region in the xy plane: closed loops, each a chain of curves running head to tail (z is
/// ignored). Loops are nested: a loop inside an even number of others is an outer boundary, one
/// inside an odd number a hole. Input loops may run either way; results are oriented with the
/// region on the left (outer loops counter-clockwise, holes clockwise).
using Region = std::vector<std::vector<NurbsCurve>>;

/// Union, intersection or difference (a minus b) of two planar regions. The curves stay exact:
/// they are split at their intersections and the pieces bounding the result are kept as they
/// are (shared boundary pieces of touching regions cancel). `tolerance` is the model-space
/// tolerance for intersections and for chaining. An empty result is an empty region.
[[nodiscard]] Region region_boolean(const Region& a, const Region& b, BooleanOp op,
                                    double tolerance = 1e-9);

/// Area of a region (outer loops minus holes), by Green's theorem on the exact curves.
[[nodiscard]] double region_area(const Region& region);

/// The region oriented with the material on the left: outer loops counter-clockwise, holes
/// clockwise (by nesting depth).
[[nodiscard]] Region oriented(const Region& region);

}  // namespace nf
