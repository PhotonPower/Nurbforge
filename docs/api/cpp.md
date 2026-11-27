# C++ API overview

Include everything with `#include <nurbforge/nurbforge.hpp>`; namespace `nf`.

```cpp
#include <nurbforge/nurbforge.hpp>
#include <numbers>

auto meridian = nf::make_arc({0,0,0}, {0,0,-1}, {1,0,0}, 10.0, 0.0, std::numbers::pi);
auto sphere   = nf::revolve(meridian, {0,0,0}, {0,0,1}, 2*std::numbers::pi);
nf::Vec3 p    = sphere.evaluate(0.25, 0.5);
```

Use from another CMake project:

```cmake
find_package(Nurbforge REQUIRED)          # after `cmake --install`
target_link_libraries(my_solver PRIVATE Nurbforge::core)
```

or with `add_subdirectory(Nurbforge)` / `FetchContent` and `-DNF_BUILD_PYTHON=OFF`.

## Solids and mass properties

`topology/solids.hpp` builds closed shells with exact NURBS faces and outward normals;
`Shell::mass_properties()` integrates volume, centroid and inertia tensor over the faces
(divergence theorem, voids subtracted), and `Shell::transformed()` maps a shell by any affine
transform:

```cpp
#include <nurbforge/topology/boolean.hpp>
#include <nurbforge/topology/solids.hpp>
nf::Shell part = nf::boolean(nf::make_box({-2, -2, 0}, {4, 4, 1}),
                             nf::make_cylinder({0, 0, -1}, {0, 0, 1}, 1.0, 3.0),
                             nf::BooleanOp::Difference);
nf::MassProperties mp = part.mass_properties();  // volume 16 - pi, centroid (0, 0, 0.5)
// mp.inertia: row-major 3x3 about the centroid (unit density); mp.principal_moments / _axes
nf::Shell turned = part.transformed(nf::Transform::rotation({0, 0, 0}, {1, 0, 0}, 0.5));
```

## Slicing shells

`analysis/shell_section.hpp` cuts closed B-rep shells into planar sections, e.g. to slice a body into
layers along a stacking direction:

```cpp
#include <nurbforge/analysis/shell_section.hpp>
const nf::Vec3 dir{0, 0, 1};
for (double h : nf::critical_heights(shell, dir)) { /* topology changes only here */ }
bool prism = nf::is_extrusion_band(shell, dir, h0, h1);       // all sections in (h0, h1) equal
nf::PlanarSection s = nf::section(shell, {{0, 0, 0.5 * (h0 + h1)}, dir});
for (const auto& loop : s.loops)            // solid on the left: area > 0 outer, < 0 hole
  for (const auto& c : loop.curves)         // c.curve in (s, t, 0), c.normal(t) outward (3D)
    ...;
```

Faces whose height along the normal depends on one parameter (extrusions and revolutions about the
normal, planes) give exact isoparametric section curves (`SectionCurve::exact`); others are traced
and fitted to `SectionOptions::fit_tolerance` (default 1e-9 of the model size). Planes through
critical heights (a plane containing a face) are rejected. All these queries are `const` and
thread-safe (checked under ThreadSanitizer in CI).

### Hatching a layer

`analysis/hatch.hpp` fills the region bounded by section loops with writing paths – parallel hatch
segments and inward contour lines, as written layer by layer in two-photon lithography:

```cpp
#include <nurbforge/analysis/hatch.hpp>
nf::HatchOptions o;
o.phase = 0.0;                 // lines n . p = phase + k * spacing, n = (-d.y, d.x)
o.contours = {0.1, 0.2};       // inward offsets of every loop
nf::Hatching h = nf::hatch(s.loops, {1, 0, 0}, 0.1, o);
for (const auto& seg : h.segments)   // sorted by line, then along the direction
  ...;                               // seg.start, seg.end, seg.line, seg.region (outer loop)
for (const auto& c : h.contours)     // c.curves: closed chain, region on the left
  ...;
```

The lines are intersected exactly with the loop curves (sign changes of the rational height of each
Bézier piece); inside and outside follow the loop orientation (winding number), so holes and
islands in holes are handled. Lines that only touch a loop or run along an edge without crossing it
are not split, so tangent lines give no zero-length segments. Contours are the boundaries of the region
shrunk by each distance: the raw offsets (`offset_curve` of the loop pieces, arcs about reflex
corners, cusp parts left out) are split at all their intersections and the fragments at the full
distance from the loops are chained. So contours are trimmed at convex corners, split at necks
narrower than twice the distance, merge where holes come close (`HatchContour::loops` lists the
source loops), and vanish where the region is too thin.

## Silhouettes, hidden lines and outlines

`analysis/silhouette.hpp` and `analysis/drawing.hpp` compute what a viewer sees of a model, for
2D drawings and figures:

```cpp
#include <nurbforge/analysis/drawing.hpp>

const nf::View view = nf::View::parallel({0, 1, 0});            // or View::perspective(eye)
for (const nf::SilhouetteCurve& c : nf::silhouette(shell, view))   // per face: c.face
  use(c.curve, c.pcurve, c.cusps);

// hidden-line removal for any curves against occluding shells
for (const nf::CurvePiece& p : nf::visibility(curves, view, {shell}))
  (p.visible ? solid : dashed).push_back(p.curve);

nf::DrawingOptions opt;
opt.s_axis = {1, 0, 0};                                            // x to the right
const nf::Drawing d = nf::make2d({shell}, {0, 1, 0}, opt);         // d.visible, d.hidden: (s, t, 0)
const nf::Outline o = nf::outline({shell}, {0, 0, -1}, opt);       // o.loops: a nf::Region
```

Silhouettes are the zero set of N · w (w the unit vector to the viewer): sign changes are
searched on a grid over the parameter domain (`SilhouetteOptions::resolution` cells per knot
span), and the traced points are moved onto the exact zero set and fitted to the tolerance.
A point of a curve is hidden if the ray from it to the viewer crosses a face, i.e. if the point
lies more than the tolerance off the tangent plane of the face where the ray meets it; curves
lying on a face are therefore not hidden by it. `make2d` and `outline` are for parallel views:
the projection onto the drawing plane is then an exact affine map of the control points.

## Shrink wrap (triangle meshes)

`mesh/shrink_wrap.hpp` builds one closed triangle mesh around faces that overlap, leave gaps or
are open. It is the only part of the kernel whose result is approximate (accurate to about the
resolution) rather than exact geometry; it uses the standard library only.

```cpp
#include <nurbforge/mesh/shrink_wrap.hpp>

nf::ShrinkWrapOptions opt;
opt.offset = 0.5;                         // distance of the hull from the input (default 0)
const nf::TriangleMesh hull = nf::shrink_wrap(std::vector<nf::Shell>{a, b}, 0.2, opt);
hull.is_watertight(); hull.volume(); hull.min_quality();   // nodes, triangles: plain vectors
nf::write_stl(hull, "envelope.stl");
```

There is an overload for loose faces (`std::vector<TrimmedSurface>`). `TriangleMesh`
(`mesh/triangle_mesh.hpp`) holds `nodes` and `triangles` and measures itself (`area`, `volume`,
`quality`, `is_watertight`, `is_consistently_oriented`); `improve_quality(mesh, min_quality)`
is the smoothing `shrink_wrap` applies when the worst triangle is below
`ShrinkWrapOptions::min_quality`. Behaviour and limits are described in
[meshing.md](../meshing.md#shrink-wrap-a-closed-hull-around-anything); remeshing with gmsh is
available from Python only.

Reference documentation: `doxygen docs/Doxyfile` → `build/docs/html/index.html`.

The same Doxygen output (XML) is rendered on the documentation site: [C++ API reference](cpp-reference.md).
