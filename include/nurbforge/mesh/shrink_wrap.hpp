#pragma once
/// @file shrink_wrap.hpp
/// Shrink wrap: one closed triangle mesh around arbitrary input (approximate, mesh level).

#include <cstddef>
#include <vector>

#include "nurbforge/mesh/triangle_mesh.hpp"
#include "nurbforge/topology/shell.hpp"
#include "nurbforge/topology/trimmed_surface.hpp"

namespace nf {

struct ShrinkWrapOptions {
  double offset{0.0};                 ///< distance of the envelope from the input (>= 0)
  bool fill_holes{true};              ///< fill everything that cannot be reached from outside
  std::size_t max_nodes{20'000'000};  ///< limit for the size of the grid
  /// Smooth the mesh until its worst triangle has this quality (TriangleMesh::quality,
  /// improve_quality); 0: return the triangles as cut out of the grid.
  double min_quality{0.3};
};

/// One closed, outward-oriented triangle mesh around faces that may overlap, leave gaps or be
/// open. Unlike everything else in the kernel the result is not exact geometry: it is accurate
/// to about `resolution`, the cell size of the grid it is built on and about the edge length of
/// the mesh. Openings and gaps narrower than about 1.5 x resolution (plus twice the offset) are
/// closed, smaller details are smoothed away.
///
/// Method: the faces are sampled into points a quarter of the resolution apart; the unsigned
/// distance to them is computed on a regular grid in a narrow band; the nodes that can be
/// reached from the border of the grid without coming closer to the input than the envelope
/// distance are outside; the surface between them and the rest is extracted with marching
/// tetrahedra, which always gives a closed, consistently oriented 2-manifold. That envelope
/// is watertight on the grid only at a distance of at least 0.75 x resolution; for a smaller
/// offset the surface that lies the difference further in is taken, from a second distance
/// field measured from the envelope (a morphological closing). Parts thinner than about
/// 1.5 x resolution do not survive this: an open surface needs an offset of at least
/// 0.75 x resolution.
///
/// With fill_holes the mesh is the outer envelope only; without it, enclosed voids larger than
/// the envelope distance are kept as inner components with normals pointing into the void
/// (the inside of a closed shell is such a void: the input is a set of surfaces). Parts
/// further apart than the gap size give separate closed components.
///
/// Throws GeometryError if resolution <= 0, offset < 0, there is nothing to wrap, the grid
/// would exceed max_nodes, or nothing is left of the envelope.
[[nodiscard]] TriangleMesh shrink_wrap(const std::vector<TrimmedSurface>& faces, double resolution,
                                       const ShrinkWrapOptions& options = {});

/// Shrink wrap around the faces of shells.
[[nodiscard]] TriangleMesh shrink_wrap(const std::vector<Shell>& shells, double resolution,
                                       const ShrinkWrapOptions& options = {});

}  // namespace nf
