add_executable(nf_bench bench.cpp)
target_link_libraries(nf_bench PRIVATE Nurbforge::core)
nf_set_warnings(nf_bench)

add_executable(nf_stress stress.cpp)
target_link_libraries(nf_stress PRIVATE Nurbforge::core)
nf_set_warnings(nf_stress)

if(NF_BUILD_TESTS)
  # smoke test: one repetition of every case (catches crashes and regressions to "hang")
  add_test(NAME benchmark_smoke COMMAND nf_bench --quick)
  # Boolean robustness: one random case per surface family (full runs: nf_stress 20)
  add_test(NAME stress_smoke COMMAND nf_stress 1)
endif()
