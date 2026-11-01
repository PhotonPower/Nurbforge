cmake_minimum_required(VERSION 3.24)

project(Nurbforge
  VERSION 0.1.0
  DESCRIPTION "Headless NURBS-based CAD kernel for scientific software"
  HOMEPAGE_URL "https://github.com/PhotonPower/Nurbforge"
  LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option(NF_BUILD_TESTS            "Build C++ unit tests"            ON)
option(NF_BUILD_PYTHON           "Build Python bindings"           ON)
option(NF_BUILD_EXAMPLES         "Build C++ examples"              ON)
option(NF_WARNINGS_AS_ERRORS     "Treat compiler warnings as errors" OFF)

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
  src/modeling/sweep.cpp
  src/modeling/gordon.cpp
  src/analysis/closest_point.cpp
  src/analysis/properties.cpp
  src/fitting/fit.cpp
  src/optics/surfaces.cpp)
add_library(Nurbforge::core ALIAS nurbforge_core)
set_target_properties(nurbforge_core PROPERTIES EXPORT_NAME core)

target_include_directories(nurbforge_core PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/generated>
  $<INSTALL_INTERFACE:include>)
nf_set_warnings(nurbforge_core)

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

# ---------------------------------------------------------------- install
if(NOT SKBUILD)
  include(GNUInstallDirs)
  install(TARGETS nurbforge_core EXPORT NurbforgeTargets
          ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
  install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
          FILES_MATCHING PATTERN "*.hpp")
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/generated/nurbforge/version.hpp
          DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/nurbforge)
  install(EXPORT NurbforgeTargets FILE NurbforgeConfig.cmake NAMESPACE Nurbforge::
          DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Nurbforge)
endif()
