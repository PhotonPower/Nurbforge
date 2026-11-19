#pragma once
// Internal: geodesics on surfaces (not installed).

#include <array>

#include "nurbforge/core/vec.hpp"
#include "nurbforge/nurbs/surface.hpp"

namespace nf::detail {

/// End point (u, v) of the geodesic of the given length that starts at (u, v) in the 3D
/// direction `dir` (projected into the tangent plane): classical Runge-Kutta on the geodesic
/// equations, steps doubled until the end point moves by less than 0.05 * tolerance. `steps`
/// is the starting step count and receives the converged one (neighbouring geodesics need
/// similar counts). Throws GeometryError if the geodesic leaves the parameter domain or the
/// surface is degenerate on the way.
[[nodiscard]] std::array<double, 2> geodesic_point(const NurbsSurface& S, double u, double v,
                                                   const Vec3& dir, double length, double tolerance,
                                                   int& steps);

}  // namespace nf::detail
