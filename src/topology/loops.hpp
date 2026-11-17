#pragma once
// Internal helpers shared by the topology modules (not installed).

#include <vector>

#include "nurbforge/nurbs/curve.hpp"
#include "nurbforge/topology/trimmed_surface.hpp"

namespace nf::detail {

/// The region bounded by curves (z = 0) with the region on their left – outer loops
/// counter-clockwise, holes clockwise – all in one parameter plane: T-junctions are resolved,
/// pieces running along each other in opposite directions cancel, the rest is chained into
/// closed loops (consecutive ends made identical). Empty if the curves do not form closed loops.
[[nodiscard]] std::vector<TrimLoop> union_loops(std::vector<NurbsCurve> curves, double eps);

}  // namespace nf::detail
