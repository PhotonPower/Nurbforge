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

Reference documentation: `doxygen docs/Doxyfile` → `build/docs/html/index.html`.
