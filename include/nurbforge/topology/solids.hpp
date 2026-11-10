#pragma once

/// @file solids.hpp
/// Closed solid primitives as B-rep shells with exact NURBS faces and outward normals: box,
/// sphere, cylinder, cone (frustum) and torus. Combine them with boolean(), map them with
/// Shell::transformed(), query them with Shell::volume() / mass_properties() / point_in_shell().

#include "nurbforge/core/vec.hpp"
#include "nurbforge/topology/shell.hpp"

namespace nf {

/// Axis-aligned box with one corner at `corner` and edge lengths size.x, size.y, size.z (> 0):
/// six planar faces.
[[nodiscard]] Shell make_box(const Vec3& corner, const Vec3& size);

/// Sphere: one revolved face (semicircle about the axis through the centre parallel to z), poles
/// at centre +- radius z.
[[nodiscard]] Shell make_sphere(const Vec3& center, double radius);

/// Circular cylinder with the base disc centred at `base`, axis direction `axis` (normalised) and
/// the given height: the extruded circle plus two planar lids trimmed by exact circles.
[[nodiscard]] Shell make_cylinder(const Vec3& base, const Vec3& axis, double radius, double height);

/// Truncated cone (frustum): radius r0 at `base`, r1 at base + height axis. One of the radii may
/// be 0 (a cone with its apex as a pole of the revolved side face; no lid there).
[[nodiscard]] Shell make_cone(const Vec3& base, const Vec3& axis, double r0, double r1,
                              double height);

/// Torus about the axis through `center` with direction `axis`: major radius R (centre to tube
/// centre) and tube radius r, 0 < r < R. One doubly periodic revolved face.
[[nodiscard]] Shell make_torus(const Vec3& center, const Vec3& axis, double major_radius,
                               double minor_radius);

}  // namespace nf
