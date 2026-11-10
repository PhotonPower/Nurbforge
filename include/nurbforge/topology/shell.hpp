#pragma once
/// @file shell.hpp
/// B-rep shells: trimmed faces sewn along shared edges.

#include <array>
#include <cstddef>
#include <string>
#include <vector>

#include "nurbforge/geometry/transform.hpp"
#include "nurbforge/nurbs/curve.hpp"
#include "nurbforge/topology/trimmed_surface.hpp"

namespace nf {

/// Use of an edge by a face boundary loop.
struct Coedge {
  std::size_t face{0};
  std::size_t loop{0};      ///< 0: outer loop, k > 0: hole k-1
  std::size_t index{0};     ///< position of the trim curve in the loop
  bool reversed{false};     ///< loop direction opposite to the edge curve
};

struct ShellEdge {
  NurbsCurve curve;              ///< 3D curve (from the first face that uses it)
  std::size_t v0{0}, v1{0};      ///< start / end vertex (equal for closed edges)
  std::vector<Coedge> uses;
};

/// Volume, centroid and inertia of the solid bounded by a closed shell, for unit density
/// (multiply volume and inertia by the density for mass and mass moments of inertia).
struct MassProperties {
  double volume{0.0};
  Vec3 centroid;
  /// Inertia tensor about the centroid, row-major: I_ij = int (|r|^2 delta_ij - r_i r_j) dV with
  /// r measured from the centroid (diagonal: moments of inertia, off-diagonal: products of
  /// inertia with the minus sign included).
  std::array<double, 9> inertia{};
  /// Eigenvalues of the inertia tensor, ascending, and the corresponding unit axes.
  std::array<double, 3> principal_moments{};
  std::array<Vec3, 3> principal_axes{};
};

/// Faces sewn along shared boundary curves. Degenerate boundary pieces (poles of revolved
/// surfaces) are not edges. Each face's boundary loops are mapped to 3D and matched: two
/// boundary curves form the same edge if each lies on the other within the sewing tolerance
/// (same or opposite direction; closed curves may start at different points). Curves that
/// overlap only partially (T-junctions) are not split.
class Shell {
 public:
  static Shell sew(std::vector<TrimmedSurface> faces, double tolerance = 1e-6);

  [[nodiscard]] const std::vector<TrimmedSurface>& faces() const noexcept { return faces_; }
  [[nodiscard]] const std::vector<Vec3>& vertices() const noexcept { return vertices_; }
  [[nodiscard]] const std::vector<ShellEdge>& edges() const noexcept { return edges_; }

  /// Edges used by exactly one face loop (holes in the shell).
  [[nodiscard]] std::vector<std::size_t> free_edges() const;
  /// Edges used more than twice.
  [[nodiscard]] std::vector<std::size_t> non_manifold_edges() const;
  /// Watertight: every edge used exactly twice.
  [[nodiscard]] bool is_closed() const;
  /// Every edge used twice is traversed in opposite directions by its two uses (normals of
  /// adjacent faces point to the same side).
  [[nodiscard]] bool is_consistently_oriented() const;

  /// Connected components (faces linked through shared edges), each as its own shell – e.g.
  /// the outer boundary and the void boundaries of a hollow solid, or the separate lumps of a
  /// Boolean result.
  [[nodiscard]] std::vector<Shell> components() const;

  /// Flips faces so that the orientation becomes consistent within every connected component.
  /// Closed components are then oriented by nesting: a component lying inside an even number
  /// of other closed components is an outer boundary (normals outwards, positive volume), one
  /// inside an odd number bounds a void (normals into the void, negative volume) – so a hollow
  /// solid keeps its void and the total volume is the solid volume.
  [[nodiscard]] Shell oriented() const;

  [[nodiscard]] double area(double rel_tol = 1e-8) const;
  /// Enclosed volume by the divergence theorem, V = 1/3 sum of the face integrals of
  /// S . (S_u x S_v). Requires a closed, consistently oriented shell (throws otherwise);
  /// positive for outward normals.
  [[nodiscard]] double volume(double rel_tol = 1e-8) const;
  /// Volume, centroid and inertia tensor by the divergence theorem (face integrals of
  /// polynomial fields in the coordinates relative to the centre of the bounding box, combined
  /// with the parallel-axis theorem). Voids (inward-oriented components) are subtracted. Same
  /// requirements as volume(); throws GeometryError if the volume is not positive.
  [[nodiscard]] MassProperties mass_properties(double rel_tol = 1e-8) const;

  /// The shell mapped by an affine transformation: every face surface is transformed (trim loops
  /// unchanged) and the faces are sewn again. Orientation-reversing maps (mirrors, negative
  /// scale factors) flip the faces, so outward normals stay outward.
  [[nodiscard]] Shell transformed(const Transform& t) const;

  /// The same shell with face f labelled labels[f] (TrimmedSurface::label; one label per face,
  /// "" for none). Labels survive Booleans, so faces can be named before an operation and found
  /// again afterwards (boundary conditions, physical groups for meshing).
  [[nodiscard]] Shell with_face_labels(const std::vector<std::string>& labels) const;
  /// Indices of the faces with the given label.
  [[nodiscard]] std::vector<std::size_t> faces_labeled(const std::string& label) const;

 private:
  [[nodiscard]] std::vector<std::vector<std::size_t>> face_components() const;

  std::vector<TrimmedSurface> faces_;
  std::vector<Vec3> vertices_;
  std::vector<ShellEdge> edges_;
  double tol_{1e-6};
};

}  // namespace nf
