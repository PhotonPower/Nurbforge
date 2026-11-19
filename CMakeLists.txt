cmake_minimum_required(VERSION 3.24)

project(Nurbforge
  VERSION 0.15.0
  DESCRIPTION "Headless NURBS-based CAD kernel for scientific software"
  HOMEPAGE_URL "https://github.com/PhotonPower/Nurbforge"
  LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# macOS: the STEP writer uses std::to_chars for double, which Apple's libc++ provides from
# macOS 13.3; release wheels target 14.0 (pyproject.toml). Fail early instead of with an
# "unavailable" compile error deep inside src/io/step.cpp.
if(APPLE AND CMAKE_OSX_DEPLOYMENT_TARGET AND CMAKE_OSX_DEPLOYMENT_TARGET VERSION_LESS 13.3)
  message(FATAL_ERROR "Nurbforge needs macOS >= 13.3 (CMAKE_OSX_DEPLOYMENT_TARGET / "
                      "MACOSX_DEPLOYMENT_TARGET is ${CMAKE_OSX_DEPLOYMENT_TARGET})")
endif()

option(NF_BUILD_TESTS            "Build C++ unit tests"            ON)
option(NF_BUILD_PYTHON           "Build Python bindings"           ON)
option(NF_BUILD_EXAMPLES         "Build C++ examples"              ON)
option(NF_WARNINGS_AS_ERRORS     "Treat compiler warnings as errors" OFF)
option(NF_BUILD_BENCHMARKS      "Build the benchmark executable"  OFF)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(CompilerWarnings)

configure_file(include/nurbforge/version.hpp.in
               ${CMAKE_CURRENT_BINARY_DIR}/generated/nurbforge/version.hpp @ONLY)

# ---------------------------------------------------------------- core library
add_library(nurbforge_core STATIC
  src/core/vec.cpp
  src/nurbs/knot_vector.cpp
  src/nurbs/basis.cpp
  src/nurbs/curve.cpp
  src/nurbs/surface.cpp
  src/nurbs/refine.cpp
  src/nurbs/reduce.cpp
  src/nurbs/join.cpp
  src/geometry/transform.cpp
  src/geometry/primitives.cpp
  src/modeling/construct.cpp
  src/modeling/pattern.cpp
  src/modeling/extend.cpp
  src/modeling/loft.cpp
  src/modeling/surface_tools.cpp
  src/modeling/sweep.cpp
  src/modeling/gordon.cpp
  src/modeling/offset.cpp
  src/modeling/fillet.cpp
  src/modeling/blend.cpp
  src/modeling/project.cpp
  src/analysis/closest_point.cpp
  src/analysis/properties.cpp
  src/analysis/intersect.cpp
  src/analysis/surface_intersect.cpp
  src/analysis/section.cpp
  src/analysis/shell_section.cpp
  src/analysis/hatch.cpp
  src/fitting/fit.cpp
  src/optics/surfaces.cpp
  src/topology/trimmed_surface.cpp
  src/topology/shell.cpp
  src/topology/arrangement.cpp
  src/topology/boolean.cpp
  src/topology/edge_blend.cpp
  src/topology/solids.cpp
  src/topology/solid_ops.cpp
  src/topology/region.cpp
  src/io/step.cpp
  src/io/step_read.cpp)
add_library(Nurbforge::core ALIAS nurbforge_core)
set_target_properties(nurbforge_core PROPERTIES EXPORT_NAME core)

target_include_directories(nurbforge_core PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/generated>
  $<INSTALL_INTERFACE:include>)
nf_set_warnings(nurbforge_core)

# std::thread (parallel seed search in the surface/surface intersection)
find_package(Threads REQUIRED)
target_link_libraries(nurbforge_core PUBLIC Threads::Threads)

if(MSVC)
  target_compile_definitions(nurbforge_core PUBLIC _USE_MATH_DEFINES NOMINMAX)
endif()

# ---------------------------------------------------------------- subprojects
if(NF_BUILD_PYTHON)
  add_subdirectory(python/bindings)
endif()

if(NF_BUILD_TESTS AND NOT SKBUILD)
  enable_testing()
  add_subdirectory(tests/cpp)
endif()

if(NF_BUILD_EXAMPLES AND NOT SKBUILD)
  add_subdirectory(examples/cpp)
endif()

if(NF_BUILD_BENCHMARKS AND NOT SKBUILD)
  add_subdirectory(benchmarks)
endif()

# ---------------------------------------------------------------- install
if(NOT SKBUILD)
  include(GNUInstallDirs)
  install(TARGETS nurbforge_core EXPORT NurbforgeTargets
          ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
  install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
          FILES_MATCHING PATTERN "*.hpp")
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/generated/nurbforge/version.hpp
          DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/nurbforge)
  # find_package(Nurbforge): config file with the dependencies, version file and targets
  include(CMakePackageConfigHelpers)
  set(NF_CMAKE_DIR ${CMAKE_INSTALL_LIBDIR}/cmake/Nurbforge)
  install(EXPORT NurbforgeTargets FILE NurbforgeTargets.cmake NAMESPACE Nurbforge::
          DESTINATION ${NF_CMAKE_DIR})
  configure_package_config_file(cmake/NurbforgeConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/NurbforgeConfig.cmake INSTALL_DESTINATION ${NF_CMAKE_DIR})
  write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/NurbforgeConfigVersion.cmake
    COMPATIBILITY SameMinorVersion)
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/NurbforgeConfig.cmake
                ${CMAKE_CURRENT_BINARY_DIR}/NurbforgeConfigVersion.cmake
          DESTINATION ${NF_CMAKE_DIR})
endif()
