#pragma once
/// @file solid_ops.hpp
/// Operations that build or modify solids for simulation: capping openings, thickening, ...

#include <string>

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

}  // namespace nf
