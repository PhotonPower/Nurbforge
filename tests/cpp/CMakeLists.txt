include(FetchContent)
find_package(Catch2 3 QUIET)
if(NOT Catch2_FOUND)
  FetchContent_Declare(Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG        v3.7.1
    GIT_SHALLOW    TRUE)
  FetchContent_MakeAvailable(Catch2)
  list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
endif()
include(Catch)

add_executable(nurbforge_tests
  test_basis.cpp
  test_curve.cpp
  test_surface.cpp
  test_transform.cpp
  test_construct.cpp
  test_refine.cpp
  test_elevate.cpp
  test_reduce.cpp
  test_split_join.cpp
  test_derivatives.cpp
  test_closest_point.cpp)
target_link_libraries(nurbforge_tests PRIVATE Nurbforge::core Catch2::Catch2WithMain)
catch_discover_tests(nurbforge_tests)
