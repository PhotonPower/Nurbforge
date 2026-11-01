#pragma once
/// @file boolean.hpp
/// Boolean operations on closed shells and point-in-solid classification.

#include "nurbforge/topology/shell.hpp"

namespace nf {

enum class BooleanOp { Union, Intersection, Difference };

/// True if p lies inside the closed shell (ray casting: intersections of rays with the trimmed
/// faces are counted; rays hitting a face boundary or grazing a face are discarded and the
/// result is the majority of several valid rays).
[[nodiscard]] bool point_in_shell(const Shell& shell, const Vec3& p, double tolerance = 1e-6);

/// Boolean operation of two closed shells. Every pair of faces is intersected; the
/// intersection curves are clipped to the part lying in both trimmed faces, the faces are split
/// along them (planar arrangement in their parameter domains), every piece is classified as
/// inside or outside the other shell, and the selected pieces are sewn:
///   Union:        A outside B + B outside A
///   Intersection: A inside B  + B inside A
///   Difference:   A outside B + (B inside A, flipped)
/// Coincident faces (detected by sampling before the surface/surface intersection) are cut along
/// each other's boundary; their common part is classified "on" with the same or the opposite
/// outward normal: kept once for union and intersection if the normals agree, kept from A for a
/// difference if they are opposite, dropped otherwise. An empty result is an empty shell.
/// The inputs are oriented outwards first. Limitations: tangential contacts, coincident faces
/// whose common part crosses the seam of a closed surface. Throws GeometryError for open input
/// shells.
[[nodiscard]] Shell boolean(const Shell& a, const Shell& b, BooleanOp op, double tolerance = 1e-6);

}  // namespace nf
