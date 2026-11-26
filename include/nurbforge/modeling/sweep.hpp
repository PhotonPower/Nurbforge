#pragma once
/// @file sweep.hpp
/// Swept surfaces: a profile curve moved along a trajectory (one rail) or between two rails.

#include <cstddef>
#include <functional>
#include <vector>

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
  double twist{0.0};       ///< additional rotation about the tangent at the end (radians),
                           ///< growing linearly with the trajectory parameter
  double tolerance{1e-6};  ///< deviation from the ideal swept surface (model units),
                           ///< estimated at test points between the placed sections
  std::size_t max_sections{1025};
};

/// Exact translational sweep S(u, v) = C(u) + T(v) - T(v0): tensor product with
/// P_ij = P_i + Q_j and w_ij = w_i w_j (works for rational curves). u: profile, v: trajectory.
[[nodiscard]] NurbsSurface sweep_translational(const NurbsCurve& profile,
                                               const NurbsCurve& trajectory);

/// Sweeps the profile along the trajectory. The profile is given in world coordinates at the
/// start of the trajectory; at parameter v it is moved to T(v) and rotated by the frame change
/// from T(v0) to T(v) (plus twist). The surface (u: profile on [0, 1], v: trajectory domain) is
/// skinned through exactly placed copies of the profile; their number is doubled until the
/// deviation from the ideal sweep is below options.tolerance. Throws GeometryError if the
/// tolerance cannot be met with max_sections or the trajectory has a vanishing tangent.
/// With SweepFrame::Fixed and no twist the exact translational sweep is returned.
[[nodiscard]] NurbsSurface sweep(const NurbsCurve& profile, const NurbsCurve& trajectory,
                                 const SweepOptions& options = {});

/// Sweep with several profiles: profile k is given in world coordinates where it crosses the
/// trajectory, at trajectory parameter params[k] (empty: the parameter of the trajectory point
/// nearest to the profile). Every profile is moved back to the trajectory start by the inverse
/// of its placement (frame as in sweep()); the moved profiles are lofted over the parameters
/// (cubic for four or more), so the section at v is that blended shape placed at v. Before the
/// first and after the last profile the shape is constant. Closed profiles should have aligned
/// seams. Same skinning and tolerance control as sweep(); u on [0, 1], v the trajectory domain.
[[nodiscard]] NurbsSurface sweep(const std::vector<NurbsCurve>& profiles,
                                 const NurbsCurve& trajectory, std::vector<double> params = {},
                                 const SweepOptions& options = {});

/// Rail revolve: the profile is revolved about the axis while its end nearer to the rail follows
/// the rail. At rail parameter v the profile is rotated about the axis to the angle of
/// rail(v) and scaled about the axis in the radial direction so that this end reaches the
/// rail's distance from the axis; with `scale_height` it is also scaled along the axis (about
/// the axis point) to the rail's height, otherwise heights are kept. The other end of a profile
/// starting on the axis stays there, e.g. a reflector with a non-round rim. u: profile on
/// [0, 1], v: rail domain (closed rail: closed surface). Same skinning and tolerance control as
/// sweep(). Throws GeometryError if the profile end or the rail touches the axis.
[[nodiscard]] NurbsSurface rail_revolve(const NurbsCurve& profile, const NurbsCurve& rail,
                                        const Vec3& axis_point, const Vec3& axis_dir,
                                        bool scale_height = false,
                                        const SweepOptions& options = {});

/// Scaling of the sections of a two-rail sweep.
enum class Sweep2Scaling {
  Uniform,     ///< the whole section is scaled by the ratio of the rail distances
  NonUniform,  ///< only the rail-to-rail direction is scaled, heights above it are kept
};

/// Two-rail sweep: the profile, given at the rail starts, is moved, rotated and uniformly
/// scaled so that the points it had on rail1(v0) and rail2(v0) follow rail1(v) and rail2(v)
/// (rail2 is linearly reparametrised to rail1's domain). The orientation about the rail-to-rail
/// axis follows the tangent of rail1. Same skinning and tolerance control as sweep().
[[nodiscard]] NurbsSurface sweep2(const NurbsCurve& profile, const NurbsCurve& rail1,
                                  const NurbsCurve& rail2, const SweepOptions& options = {});

/// Two-rail sweep with several profiles and a choice of scaling. Profile k is given in world
/// coordinates where it spans the rails, at the parameter params[k] of rail1 (empty: the
/// parameter of the point of rail1 nearest to the profile). Every profile is moved back to the
/// rail starts by the inverse of its placement; the moved profiles are lofted over the
/// parameters (cubic for four or more), and the section at v is that blended shape placed at v,
/// as in sweep() with several profiles. Before the first and after the last profile the shape
/// is constant. With Sweep2Scaling::NonUniform the section is stretched along the rail-to-rail
/// direction only (its extent normal to that direction is kept), otherwise uniformly.
[[nodiscard]] NurbsSurface sweep2(const std::vector<NurbsCurve>& profiles, const NurbsCurve& rail1,
                                  const NurbsCurve& rail2, std::vector<double> params = {},
                                  Sweep2Scaling scaling = Sweep2Scaling::Uniform,
                                  const SweepOptions& options = {});

/// Pipe surface: a circle swept along `path` with the rotation-minimising frame, its radius
/// given by `radius(t)` for the normalised path parameter t in [0, 1] (t = 0 at the path start).
/// u runs around the circle on [0, 1], v along the path domain. Same skinning and tolerance
/// control as sweep(); the radius function must be smooth and positive.
[[nodiscard]] NurbsSurface pipe_surface(const NurbsCurve& path,
                                        const std::function<double(double)>& radius,
                                        const SweepOptions& options = {});

}  // namespace nf
