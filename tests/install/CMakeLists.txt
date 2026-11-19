# Consumer project for the installed CMake package (CI: install, then configure and build this
# with CMAKE_PREFIX_PATH pointing at the install prefix). It must not need find_package(Threads).
cmake_minimum_required(VERSION 3.24)
project(NurbforgeConsumer LANGUAGES CXX)

find_package(Nurbforge 0.16 REQUIRED)

add_executable(consumer main.cpp)
target_compile_features(consumer PRIVATE cxx_std_20)
target_link_libraries(consumer PRIVATE Nurbforge::core)
