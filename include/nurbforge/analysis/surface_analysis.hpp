#pragma once
/// @file surface_analysis.hpp
/// Analysis of surfaces: curvature and draft angle maps, direction, continuity between
/// surfaces and across the edges of a shell, duplicate surfaces.

#include <cstddef>
#include <vector>

#include "nurbforge/core/vec.hpp"
#include "nurbforge/modeling/extend.hpp"
#include "nurbforge/nurbs/surface.hpp"
#include "nurbforge/topology/shell.hpp"
#include "nurbforge/topology/trimmed_surface.hpp"

namespace nf {

// ------------------------------------------------------------------ maps

/// A scalar sampled on a grid over the parameter domain of a surface.
struct SurfaceMapRange {
  double min{0.0}, max{0.0};
  double min_u{0.0}, min_v{0.0};  ///< parameters of the smallest sample
  double max_u{0.0}, max_v{0.0};  ///< parameters of the largest sample
};

/// Curvatures on a grid (curvature analysis). Values are stored row by row:
/// index = i * v.size() + j for the parameters (u[i], v[j]).
struct CurvatureMap {
  std::vector<double> u, v;
  std::vector<double> gaussian;  ///< K = k1 * k2
  std::vector<double> mean;      ///< H = (k1 + k2) / 2, sign relative to S_u x S_v
  std::vector<double> k_max;     ///< larger principal curvature
  std::vector<double> k_min;     ///< smaller principal curvature
  /// Samples inside the trimmed region (all true for an untrimmed surface); the ranges below
  /// cover these samples only.
  std::vector<bool> inside;
  SurfaceMapRange gaussian_range, mean_range, k_max_range, k_min_range;
  /// Smallest radius of curvature 1 / max(|k_max|, |k_min|) among the samples (infinity for a
  /// plane): the largest offset distance or tool radius that fits.
  double min_radius{0.0};
};

/// Curvature analysis of a surface: Gaussian, mean and principal curvatures at
/// `samples_per_span` + 1 parameters per knot span in each direction (knot lines included, so
/// extrema on them are met), with the ranges and where they are attained.
[[nodiscard]] CurvatureMap curvature_map(const NurbsSurface& surface,
                                         std::size_t samples_per_span = 8);

/// Curvature analysis of a trimmed face: the map of its surface with `inside` set from the trim.
[[nodiscard]] CurvatureMap curvature_map(const TrimmedSurface& face,
                                         std::size_t samples_per_span = 8);

/// Draft angle (radians) at a point of a surface for a pull direction: the angle between the
/// tangent plane and the direction, asin(N . d) with the unit normal N = S_u x S_v. Positive
/// where the normal has a component along the pull direction (the face can be withdrawn),
/// zero for a wall parallel to it, negative for an undercut.
[[nodiscard]] double draft_angle(const NurbsSurface& surface, double u, double v,
                                 const Vec3& pull_direction);

/// Draft angles on a grid (draft angle analysis); layout as CurvatureMap.
struct DraftAngleMap {
  std::vector<double> u, v;
  std::vector<double> angle;  ///< radians
  std::vector<bool> inside;
  SurfaceMapRange range;
};

[[nodiscard]] DraftAngleMap draft_angle_map(const NurbsSurface& surface, const Vec3& pull_direction,
                                            std::size_t samples_per_span = 8);
[[nodiscard]] DraftAngleMap draft_angle_map(const TrimmedSurface& face, const Vec3& pull_direction,
                                            std::size_t samples_per_span = 8);

// ------------------------------------------------------------------ direction

/// Parameter directions and orientation of a surface.
struct SurfaceDirection {
  Vec3 point;   ///< the surface point at the middle of the parameter domain
  Vec3 u_dir;   ///< unit direction of increasing u there
  Vec3 v_dir;   ///< unit direction of increasing v there
  Vec3 normal;  ///< unit normal S_u x S_v there
  bool closed_u{false}, closed_v{false};  ///< the sides u_min / u_max (v_min / v_max) coincide
  /// Sides that collapse to a point, in the order u_min, u_max, v_min, v_max.
  bool pole[4]{false, false, false, false};
  bool planar{false};  ///< the control points lie in a plane within the tolerance
};

/// Directions of a surface at the middle of its domain, closed directions and poles.
/// NurbsSurface::transposed() swaps u and v and so flips the normal.
[[nodiscard]] SurfaceDirection surface_direction(const NurbsSurface& surface,
                                                 double tolerance = 1e-9);

// ------------------------------------------------------------------ continuity

struct SurfaceContinuityTolerance {
  double position{1e-6};   ///< largest gap along the edge (model units)
  double angle{1e-6};      ///< largest angle between the tangent planes for G1 (radians)
  double curvature{1e-4};  ///< largest difference of the curvature tensors for G2, relative to
                           ///< the largest principal curvature (or to 1 / size of the
                           ///< surfaces where both are nearly flat)
};

/// Continuity of two surfaces along a common edge, from samples along the edge.
struct SurfaceContinuity {
  /// -1 the surfaces do not meet along the edge, 0 position (G0), 1 tangent plane (G1),
  /// 2 curvature (G2: equal normal curvatures in every direction).
  int geometric{-1};
  double max_gap{0.0};    ///< largest distance of an edge sample from the other surface
  double max_angle{0.0};  ///< largest angle between the tangent planes (radians)
  double max_curvature_difference{0.0};  ///< largest difference of the curvature tensors
};

/// Continuity of the surface `a` along its side `side` with the surface `b`: every sample of
/// the boundary curve is compared with the closest point of b (which need not be on a side of
/// b). The tangent planes are compared whatever the orientations of the two surfaces.
[[nodiscard]] SurfaceContinuity continuity(const NurbsSurface& a, SurfaceSide side,
                                           const NurbsSurface& b,
                                           const SurfaceContinuityTolerance& tolerance = {},
                                           std::size_t samples = 33);

/// Continuity across every edge of a shell: result[e] for Shell::edges()[e], from samples
/// along the edge on its two faces. Edges not used by exactly two faces (free and
/// non-manifold edges) get geometric = -1. An edge with geometric >= 1 is a smooth (tangent)
/// edge: the seam of a cylinder, the border of a fillet.
[[nodiscard]] std::vector<SurfaceContinuity> edge_continuity(
    const Shell& shell, const SurfaceContinuityTolerance& tolerance = {}, std::size_t samples = 9);

// ------------------------------------------------------------------ duplicates

/// Groups of geometrically identical surfaces (SelDup): the same point set within the
/// tolerance, whatever the parametrisation, orientation, degrees or knots. Every group lists
/// at least two indices in increasing order; the first is the one to keep. Surfaces are
/// compared by the distances of a grid of points of each to the other.
[[nodiscard]] std::vector<std::vector<std::size_t>> find_duplicate_surfaces(
    const std::vector<NurbsSurface>& surfaces, double tolerance = 1e-6);

}  // namespace nf
