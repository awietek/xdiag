set(XDIAG_EXAMPLES
  usage_examples
  user_guide
  hello_world
  spinhalf_chain_e0
  spinhalf_chain_gs_corr_symmetries
  spinhalf_chain_structure_factor
  hubbard_greens_f
  entanglement_ground_state
  tpq_shastry_sutherland
  specific_heat_randomtj
  tos_alpha_xx
  tos_triangular
)

foreach(name ${XDIAG_EXAMPLES})
  set(src "${name}/main.cpp")
  message(STATUS "Add example: ${name} ${src}")
  add_executable(${name} ${src})
  target_link_libraries(${name} PUBLIC ${XDIAG_LIBRARY})
endforeach()


# if(XDIAG_DISTRIBUTED)
#   add_executable(tests_distributed tests_distributed.cpp ${XDIAG_TEST_DISTRIBUTED_SOURCES} ${XDIAG_TESTCASES_SOURCES})
#   target_link_libraries(tests_distributed PUBLIC ${XDIAG_LIBRARY})
#   add_test(NAME XdiagTestDistributed COMMAND tests_distributed)
# endif()
