#pragma once
/// @file edge_blend.hpp
/// Fillets and chamfers of the edges of a shell (solid).

#include <cstddef>
#include <string>
#include <vector>

#include "nurbforge/core/vec.hpp"
#include "nurbforge/topology/shell.hpp"

namespace nf {

/// Fillet radius or chamfer distance at a point on the selected edges (variable size).
struct BlendSize {
  Vec3 point;
  double size{0.0};
};

/// Options of fillet_edges and chamfer_edges.
struct EdgeBlendOptions {
  /// Label of the new faces ("" : "fillet" for fillets, "chamfer" for chamfers).
  std::string label{};
  /// Sewing tolerance of the result (0: the shell's tolerance). The blend surfaces, contact
  /// curves and end curves are computed to a twentieth of it.
  double tolerance{0.0};
  /// Variable size: radius (fillets) or distance (chamfers) at points on the selected edges
  /// (within the sewing tolerance). Along each edge chain the size passes through the values of
  /// its points in chain order, blending smoothly between neighbours (quintic smooth step, zero
  /// slope and curvature at the points: s(t) = s0 + (s1 - s0)(10t^3 - 15t^4 + 6t^5), t the
  /// arc-length fraction between them); before the first and after the
  /// last point it is constant (closed chains: periodic). A chain without points uses the
  /// size argument. Where blends meet at a corner of three edges their fillet radii must agree.
  std::vector<BlendSize> sizes{};
  /// fillet_edges only: curvature-continuous (G2) blends instead of rolling-ball fillets. The
  /// blend runs between the same contact curves as the fillet of the given radius, but its
  /// cross-sections are quintic curves that leave the two faces with their tangent and with
  /// their normal curvature in the cross direction, instead of circular arcs (which jump from
  /// the curvature of the face to 1 / radius). Reflections run smoothly across such a blend.
  /// Supported for chains and closed loops whose ends meet no other blended chain.
  bool curvature_continuous{false};
};

/// Indices of the edges shared by a face labelled `label_a` and a face labelled `label_b`
/// (e.g. edges_between(box, "top", "front")).
[[nodiscard]] std::vector<std::size_t> edges_between(const Shell& shell, const std::string& label_a,
                                                     const std::string& label_b);

/// Rounds edges of a shell with a rolling-ball fillet (constant radius, or variable with
/// options.sizes).
///
/// The selected edges are grouped by their two adjacent faces; each group (one edge, a chain of
/// edges, or a closed loop such as the rim of a cylinder) is filleted with one surface: the
/// rolling-ball fillet between the two face surfaces (fillet_surfaces; the ball lies inside
/// the material at a convex edge and outside at a concave edge, which adds material). The two
/// faces are cut back along the contact curves. An open chain ends at a vertex of three faces
/// (the two faces and a third one, E):
///   - only this chain ends there: the fillet is trimmed by E along the exact intersection curve
///     and E's boundary is rerouted along it (convex: the corner is cut off; concave: E grows
///     into the corner);
///   - a second chain (between one of the two faces and E) ends there and E is tangent to the
///     other face (a tangent-continuous chain, e.g. the rim of a rounded rectangle): both
///     fillets end at their common cross-section;
///   - a second chain ends there at a sharp edge: the two fillets are trimmed by each other
///     along their exact intersection curve (mitre); where one of them runs past the edge
///     between the other faces, it is trimmed by the face there, which grows into the corner;
///   - three chains end there (all three edges of the corner): every fillet ends at its
///     cross-section through the ball touching all three faces, and the corner is closed by
///     the spherical triangle of that ball between the three contact points (exact, tangent to
///     the three fillets).
/// The new faces are labelled options.label ("fillet"); all other labels are kept.
///
/// Throws GeometryError if an edge is not shared by exactly two faces, if a chain ends at a
/// vertex that does not join exactly three faces and three edges, if a chain continues
/// tangentially into an edge that is not selected, if a convex and a concave chain meet, if the
/// radius is too large for the faces (a contact curve leaves a face through another edge), or
/// if the result is not closed.
[[nodiscard]] Shell fillet_edges(const Shell& shell, const std::vector<std::size_t>& edges,
                                 double radius, const EdgeBlendOptions& options = {});

/// Chamfers edges of a shell: a ruled cut from the curve at (geodesic) distance `distance` from
/// the edge on one face to the curve at the same distance on the other face (variable with
/// options.sizes). Grouping, ends, labels ("chamfer") and errors as for fillet_edges; convex
/// edges lose material, concave edges gain it. At a corner of three chamfered edges the
/// chamfers end at the straight lines between the points where their contact curves cross,
/// and the corner is closed by the planar triangle of these points (faces at the corner must be
/// planar).
[[nodiscard]] Shell chamfer_edges(const Shell& shell, const std::vector<std::size_t>& edges,
                                  double distance, const EdgeBlendOptions& options = {});

}  // namespace nf
