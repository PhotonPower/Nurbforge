add_executable(nf_bench bench.cpp)
target_link_libraries(nf_bench PRIVATE Nurbforge::core)
nf_set_warnings(nf_bench)

if(NF_BUILD_TESTS)
  # smoke test: one repetition of every case (catches crashes and regressions to "hang")
  add_test(NAME benchmark_smoke COMMAND nf_bench --quick)
endif()
