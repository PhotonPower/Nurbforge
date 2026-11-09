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

Reference documentation: `doxygen docs/Doxyfile` → `build/docs/html/index.html`.

The same Doxygen output (XML) is rendered on the documentation site: [C++ API reference](cpp-reference.md).
