#pragma once
/// @file extend.hpp
/// Extrapolation: extending curves and surfaces beyond their boundary.

#include "nurbforge/nurbs/curve.hpp"
#include "nurbforge/nurbs/surface.hpp"

namespace nf {

/// How the geometry is continued.
enum class Extension {
  Natural,  ///< polynomial (rational) continuation of the boundary Bézier piece: C-infinity at
            ///< the old boundary; conics stay conics. Fails if weights would become <= 0.
  Linear,   ///< straight continuation along the boundary tangent: G1, and C1 for non-rational
            ///< geometry (the parameter speed is matched).
};

enum class CurveEnd { Start, End };
enum class SurfaceSide { UMin, UMax, VMin, VMax };

/// Extends the curve at one end by the arc length `length` (> 0). The original parameters are
/// kept: the domain grows to [a - delta, b] or [a, b + delta]. The result is clamped.
/// Throws GeometryError for length <= 0 or if a natural extension of a rational curve would
/// need non-positive weights (use a shorter length or Extension::Linear).
[[nodiscard]] NurbsCurve extend(const NurbsCurve& c, CurveEnd end, double length,
                                Extension type = Extension::Natural);

/// Extends the surface across one side. `length` is the arc length of the extension measured
/// along the iso-curve through the middle of the side (other iso-curves are extended by the
/// same parameter amount). Original parameters are kept; the domain grows on that side.
/// Linear extension requires constant weights along the extension direction (non-rational
/// iso-curves across the side), e.g. a cylinder can be extended linearly along its axis.
[[nodiscard]] NurbsSurface extend(const NurbsSurface& s, SurfaceSide side, double length,
                                  Extension type = Extension::Natural);

}  // namespace nf
