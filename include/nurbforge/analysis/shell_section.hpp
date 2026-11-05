#pragma once

/// @file shell_section.hpp
/// Planar sections of B-rep shells (closed loops with orientation) and the heights at which the
/// topology of the sections changes; used to slice bodies into layers (e.g. FourierGrating's CAD
/// slicer). All functions are const queries and thread-safe (no shared mutable state).

#include <cstddef>
#include <vector>

#include "nurbforge/analysis/section.hpp"
#include "nurbforge/topology/shell.hpp"

namespace nf {

struct SectionOptions {
  /// Chaining tolerance (model units): ends of neighbouring section pieces closer than this are
  /// joined (the shell's sewing tolerance; trim loops are reconstructed to a tenth of it).
  double join_tolerance{1e-6};
  /// Accuracy of the traced and fitted section curves on general faces (model units); 0: 1e-9 of
  /// the model size. Faces whose height along the plane normal depends on one parameter only
  /// (extrusions and revolutions about the normal, planar faces) give exact isoparametric curves.
  double fit_tolerance{0.0};
  /// First in-plane axis s (projected into the plane and normalised); zero: chosen automatically.
  /// The second axis is t = normal x s.
  Vec3 s_axis{0, 0, 0};
};

/// One piece of a section loop: a curve on one face of the shell.
struct SectionCurve {
  NurbsCurve curve;      ///< in the plane coordinates (s, t, 0)
  NurbsCurve curve3d;    ///< the same curve in 3D (same parameter)
  NurbsCurve pcurve;     ///< (u, v, 0) on the surface of the source face (same parameter)
  std::size_t face{0};   ///< index of the source face in the shell
  bool exact{false};     ///< isoparametric curve of the face; otherwise traced and fitted
  NurbsSurface surface;  ///< surface of the source face (for the normal)
  double normal_sign{1.0};

  /// Outward unit normal of the shell at parameter t of the curve (surface normal of the source
  /// face, oriented out of the solid).
  [[nodiscard]] Vec3 normal(double t) const;
};

/// Closed loop of section pieces (end of one = start of the next within the join tolerance),
/// oriented with the solid on the left in (s, t): outer loops counter-clockwise (area > 0), voids
/// and holes clockwise (area < 0).
struct SectionLoop {
  std::vector<SectionCurve> curves;
  double area{0.0};  ///< signed area in (s, t)
};

struct PlanarSection {
  Vec3 origin, s_axis, t_axis, normal;  ///< plane coordinates: p = origin + s s_axis + t t_axis
  std::vector<SectionLoop> loops;
};

/// Section of a closed shell with a plane: the intersection curves of all faces, clipped to the
/// trimmed regions, chained into closed loops and oriented by the outward normals (the shell is
/// oriented first if necessary). Throws GeometryError if a chain does not close (open shell or
/// tolerance problem; the message gives the position) or if the plane contains a face (planes
/// through critical heights, see critical_heights, are not supported).
[[nodiscard]] PlanarSection section(const Shell& shell, const Plane& plane,
                                    const SectionOptions& options = {});

/// Heights h = p . dir (dir normalised) at which the topology of the sections of the shell with the
/// planes p . dir = h can change, sorted and merged within tol.length: vertices, planar faces
/// normal to dir, interior extrema of the height on faces (S_u . dir = S_v . dir = 0, Newton from
/// Bezier subdivision) and on edges (C' . dir = 0). Between two neighbouring critical heights the
/// sections have the same number of loops and pieces.
[[nodiscard]] std::vector<double> critical_heights(const Shell& shell, const Vec3& dir,
                                                   const Tolerance& tol = {});

/// True if every face of the shell that meets the open band h0 < p . dir < h1 has a normal
/// perpendicular to dir there (within tol.angular), i.e. the band is an extrusion along dir and all
/// its sections are equal. Checked at samples of the face inside the band and on its Bezier
/// control nets.
[[nodiscard]] bool is_extrusion_band(const Shell& shell, const Vec3& dir, double h0, double h1,
                                     const Tolerance& tol = {});

}  // namespace nf
