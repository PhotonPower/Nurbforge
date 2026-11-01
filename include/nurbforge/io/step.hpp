#pragma once
/// @file step.hpp
/// Native STEP (ISO 10303-21) writer for B-rep shells, AP214 (AUTOMOTIVE_DESIGN) schema.

#include <string>
#include <vector>

#include "nurbforge/topology/shell.hpp"

namespace nf {

/// Length unit declared in the STEP file. Coordinates are written unchanged: the unit says what
/// one model unit means (default: millimetre, the usual CAD convention).
enum class StepLengthUnit { Millimetre, Centimetre, Metre, Inch };

struct StepOptions {
  StepLengthUnit unit{StepLengthUnit::Millimetre};
  /// Product name (PRODUCT entity) and FILE_NAME name.
  std::string name{"nurbforge"};
  std::string author{};
  std::string organization{};
  /// ISO 8601 time stamp for FILE_NAME; empty: current UTC time.
  std::string timestamp{};
  /// Modelling uncertainty (UNCERTAINTY_MEASURE_WITH_UNIT, model units); readers such as
  /// OpenCASCADE use it as the initial vertex / edge tolerance.
  double uncertainty{1e-6};
};

/// The STEP file contents for the given shells (see write_step).
[[nodiscard]] std::string to_step(const std::vector<Shell>& shells, const StepOptions& options = {});

/// Writes the shells to an ISO 10303-21 file with the AP214 schema (AUTOMOTIVE_DESIGN), one
/// product whose shape is:
///  - closed shells: MANIFOLD_SOLID_BREP (BREP_WITH_VOIDS if a shell has inner void components)
///    in an ADVANCED_BREP_SHAPE_REPRESENTATION; shells that are not consistently oriented are
///    oriented first, and outer components get outward normals;
///  - open shells: SHELL_BASED_SURFACE_MODEL(OPEN_SHELL) in a
///    MANIFOLD_SURFACE_SHAPE_REPRESENTATION (if both kinds are present, solids and surface model
///    go into one plain SHAPE_REPRESENTATION).
/// Faces whose normal must point against the surface normal (inside-out components, void
/// shells of BREP_WITH_VOIDS) are written with same_sense .F. and bound orientation .F.
/// Faces are ADVANCED_FACEs on B_SPLINE_SURFACE_WITH_KNOTS (rational: the complex entity with
/// RATIONAL_B_SPLINE_SURFACE), bounded by EDGE_LOOPs of ORIENTED_EDGEs; every shell edge is one
/// EDGE_CURVE (its 3D B-spline curve, as stored in the shell) shared by its uses. Empty shells
/// are skipped. Throws std::runtime_error if the file cannot be written.
void write_step(const std::string& path, const std::vector<Shell>& shells, const StepOptions& options = {});

}  // namespace nf
