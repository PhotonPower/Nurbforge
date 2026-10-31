find_package(Python 3.10 REQUIRED COMPONENTS Interpreter Development.Module)

# pybind11: prefer an installed package (pip/conda/dnf), fall back to FetchContent.
if(NOT pybind11_FOUND)
  execute_process(COMMAND "${Python_EXECUTABLE}" -m pybind11 --cmakedir
                  OUTPUT_VARIABLE _pybind11_dir OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
  find_package(pybind11 2.12 CONFIG QUIET HINTS "${_pybind11_dir}")
endif()
if(NOT pybind11_FOUND)
  include(FetchContent)
  FetchContent_Declare(pybind11
    GIT_REPOSITORY https://github.com/pybind/pybind11.git
    GIT_TAG        v2.13.6
    GIT_SHALLOW    TRUE)
  FetchContent_MakeAvailable(pybind11)
endif()

pybind11_add_module(_core MODULE module.cpp)
target_link_libraries(_core PRIVATE Nurbforge::core)
target_compile_definitions(_core PRIVATE NF_VERSION_INFO="${PROJECT_VERSION}")

if(SKBUILD)
  install(TARGETS _core LIBRARY DESTINATION nurbforge)
else()
  # In-tree development: put the extension next to the Python sources.
  set_target_properties(_core PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/python/nurbforge"
    LIBRARY_OUTPUT_DIRECTORY_RELEASE "${PROJECT_SOURCE_DIR}/python/nurbforge"
    LIBRARY_OUTPUT_DIRECTORY_DEBUG "${PROJECT_SOURCE_DIR}/python/nurbforge")
endif()
