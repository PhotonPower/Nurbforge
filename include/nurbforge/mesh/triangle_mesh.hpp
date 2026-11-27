#pragma once
/// @file triangle_mesh.hpp
/// A plain triangle mesh, the result of the mesh-level tools (shrink_wrap). Everything else in
/// the kernel is exact geometry; meshes for solvers are made by gmsh from that (Python layer).

#include <array>
#include <cstddef>
#include <string>
#include <vector>

#include "nurbforge/core/vec.hpp"

namespace nf {

/// Triangle surface mesh with 0-based connectivity.
struct TriangleMesh {
  std::vector<Vec3> nodes;
  std::vector<std::array<std::size_t, 3>> triangles;

  [[nodiscard]] double area() const;
  /// Enclosed volume (divergence theorem); positive for an outward-oriented closed mesh.
  [[nodiscard]] double volume() const;
  /// Shape quality per triangle: 4 sqrt(3) area / (sum of the squared side lengths), 1 for an
  /// equilateral triangle and 0 for a degenerate one.
  [[nodiscard]] std::vector<double> quality() const;
  /// Smallest quality(); 1 for an empty mesh.
  [[nodiscard]] double min_quality() const;
  /// Every edge is shared by exactly two triangles.
  [[nodiscard]] bool is_watertight() const;
  /// Every directed edge occurs once: neighbours traverse a shared edge in opposite directions.
  [[nodiscard]] bool is_consistently_oriented() const;
};

/// Moves the nodes of a closed mesh within the surface towards the centres of their neighbours
/// until the worst triangle reaches `min_quality` or `max_iterations` steps are done (tangential
/// Laplacian smoothing; the quality is checked every five steps). The connectivity is kept, and
/// since the nodes move in their tangent planes the shape changes to second order only.
/// Returns the number of steps taken.
std::size_t improve_quality(TriangleMesh& mesh, double min_quality,
                            std::size_t max_iterations = 50);

/// Writes the mesh as a binary STL file. Throws GeometryError if the file cannot be written.
void write_stl(const TriangleMesh& mesh, const std::string& path);

}  // namespace nf
