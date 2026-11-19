#pragma once
/// @file edge_blend.hpp
/// Fillets and chamfers of the edges of a shell (solid).

#include <cstddef>
#include <string>
#include <vector>

#include "nurbforge/topology/shell.hpp"

namespace nf {

/// Options of fillet_edges and chamfer_edges.
struct EdgeBlendOptions {
  /// Label of the new faces ("" : "fillet" for fillets, "chamfer" for chamfers).
  std::string label{};
  /// Sewing tolerance of the result (0: the shell's tolerance). The blend surfaces, contact
  /// curves and end curves are computed to a twentieth of it.
  double tolerance{0.0};
};

/// Indices of the edges shared by a face labelled `label_a` and a face labelled `label_b`
/// (e.g. edges_between(box, "top", "front")).
[[nodiscard]] std::vector<std::size_t> edges_between(const Shell& shell, const std::string& label_a,
                                                     const std::string& label_b);

/// Rounds edges of a shell with a constant-radius rolling-ball fillet.
///
/// The selected edges are grouped by their two adjacent faces; each group (one edge, a chain of
/// edges, or a closed loop such as the rim of a cylinder) is filleted with one surface: the
/// rolling-ball fillet between the two face surfaces (fillet_surfaces; the ball lies inside
/// the material at a convex edge and outside at a concave edge, which adds material). The two
/// faces are cut back along the contact curves. An open chain ends at a vertex where exactly
/// one further face E meets the two faces (e.g. the side face at the end of a box edge): the
/// fillet surface is trimmed by E along the exact intersection curve, and E's boundary is
/// rerouted along that curve (convex: the corner is cut off; concave: E grows into the
/// corner). The new faces are labelled options.label ("fillet"); all other labels are kept.
///
/// Groups are processed one after another. Throws GeometryError if an edge is not shared by
/// exactly two faces, if two groups (or an edge chain and another selected edge) share a
/// vertex (vertex blends are not supported yet), if an end vertex touches more than one
/// further face, if the radius is too large for the faces (the contact curve leaves a face
/// through another edge), or if the result is not closed.
[[nodiscard]] Shell fillet_edges(const Shell& shell, const std::vector<std::size_t>& edges,
                                 double radius, const EdgeBlendOptions& options = {});

/// Chamfers edges of a shell: a ruled cut from the curve at (geodesic) distance `distance` from
/// the edge on one face to the curve at the same distance on the other face. Grouping, ends,
/// labels ("chamfer") and errors as for fillet_edges; convex edges lose material, concave edges
/// gain it.
[[nodiscard]] Shell chamfer_edges(const Shell& shell, const std::vector<std::size_t>& edges,
                                  double distance, const EdgeBlendOptions& options = {});

}  // namespace nf
