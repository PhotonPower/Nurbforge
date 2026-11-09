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

Reference documentation: `doxygen docs/Doxyfile` → `build/docs/html/index.html`.

The same Doxygen output (XML) is rendered on the documentation site: [C++ API reference](cpp-reference.md).
