#pragma once
/// @file shell.hpp
/// B-rep shells: trimmed faces sewn along shared edges.

#include <cstddef>
#include <vector>

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

  /// Flips faces so that the orientation becomes consistent (per connected component). For a
  /// closed shell the normals then point outwards (enclosed volume positive).
  [[nodiscard]] Shell oriented() const;

  [[nodiscard]] double area(double rel_tol = 1e-8) const;
  /// Enclosed volume by the divergence theorem, V = 1/3 sum of the face integrals of
  /// S . (S_u x S_v). Requires a closed, consistently oriented shell (throws otherwise);
  /// positive for outward normals.
  [[nodiscard]] double volume(double rel_tol = 1e-8) const;

 private:
  std::vector<TrimmedSurface> faces_;
  std::vector<Vec3> vertices_;
  std::vector<ShellEdge> edges_;
  double tol_{1e-6};
};

}  // namespace nf
