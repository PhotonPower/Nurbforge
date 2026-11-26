#pragma once
// Internal: clipping of parameter-space curves to trimmed regions (not installed).

#include <array>
#include <vector>

#include "nurbforge/nurbs/curve.hpp"
#include "nurbforge/topology/trimmed_surface.hpp"

namespace nf::detail {

/// Parameter intervals of the p-curve inside the trimmed region (bisection on contains()).
[[nodiscard]] std::vector<std::array<double, 2>> inside_intervals(const TrimmedSurface& f,
                                                                  const NurbsCurve& pc);

}  // namespace nf::detail
