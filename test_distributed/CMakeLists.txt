add_executable(tests
  tests.cpp
  mpi/test_dot_mpi.cpp
  basis/tj_distributed/test_basis_np.cpp
)
target_link_libraries(tests PRIVATE hydra)
add_test(NAME HydraTest COMMAND tests)
