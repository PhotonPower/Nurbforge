#pragma once
/// @file sweep.hpp
/// Swept surfaces: a profile curve moved along a trajectory (one rail) or between two rails.

#include <cstddef>

#include "nurbforge/nurbs/curve.hpp"
#include "nurbforge/nurbs/surface.hpp"

namespace nf {

/// Orientation of the profile along the trajectory.
enum class SweepFrame {
  RotationMinimizing,  ///< rotation-minimising frame (double reflection, Wang et al. 2008):
                       ///< no spurious twist, defined at inflections and straight parts
  Fixed,               ///< pure translation (no rotation)
};

struct SweepOptions {
  SweepFrame frame{SweepFrame::RotationMinimizing};
  double twist{0.0};            ///< additional rotation about the tangent at the end (radians),
                                ///< growing linearly with the trajectory parameter
  double tolerance{1e-6};       ///< deviation from the ideal swept surface (model units),
                                ///< estimated at test points between the placed sections
  std::size_t max_sections{1025};
};

/// Exact translational sweep S(u, v) = C(u) + T(v) - T(v0): tensor product with
/// P_ij = P_i + Q_j and w_ij = w_i w_j (works for rational curves). u: profile, v: trajectory.
[[nodiscard]] NurbsSurface sweep_translational(const NurbsCurve& profile, const NurbsCurve& trajectory);

/// Sweeps the profile along the trajectory. The profile is given in world coordinates at the
/// start of the trajectory; at parameter v it is moved to T(v) and rotated by the frame change
/// from T(v0) to T(v) (plus twist). The surface (u: profile on [0, 1], v: trajectory domain) is
/// skinned through exactly placed copies of the profile; their number is doubled until the
/// deviation from the ideal sweep is below options.tolerance. Throws GeometryError if the
/// tolerance cannot be met with max_sections or the trajectory has a vanishing tangent.
/// With SweepFrame::Fixed and no twist the exact translational sweep is returned.
[[nodiscard]] NurbsSurface sweep(const NurbsCurve& profile, const NurbsCurve& trajectory,
                                 const SweepOptions& options = {});

/// Two-rail sweep: the profile, given at the rail starts, is moved, rotated and uniformly
/// scaled so that the points it had on rail1(v0) and rail2(v0) follow rail1(v) and rail2(v)
/// (rail2 is linearly reparametrised to rail1's domain). The orientation about the rail-to-rail
/// axis follows the tangent of rail1. Same skinning and tolerance control as sweep().
[[nodiscard]] NurbsSurface sweep2(const NurbsCurve& profile, const NurbsCurve& rail1,
                                  const NurbsCurve& rail2, const SweepOptions& options = {});

}  // namespace nf
