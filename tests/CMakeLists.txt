set(XDIAG_TEST_SOURCES
  bits/test_bitops.cpp

  io/test_args.cpp
  io/test_file_toml.cpp
  io/test_file_h5.cpp

  combinatorics/test_binomial.cpp
  combinatorics/test_subsets.cpp
  combinatorics/test_bit_patterns.cpp
  combinatorics/test_combinations.cpp
  combinatorics/test_combinations_index.cpp
  combinatorics/test_lin_table.cpp
  combinatorics/test_fermi_table.cpp

  basis/spinhalf/test_spinhalf_basis_sublattice.cpp
  basis/spinhalf/test_spinhalf_basis.cpp
  basis/spinhalf/test_spinhalf_basis_iterator.cpp
  basis/electron/test_basis_electron.cpp
  basis/tj/test_basis_tj.cpp

  symmetries/group_action/test_group_action.cpp
  symmetries/group_action/test_group_action_lookup.cpp
  symmetries/group_action/test_group_action_sublattice.cpp
  symmetries/operations/test_symmetry_operations.cpp
  # symmetries/operations/test_representative_list_omp.cpp
  symmetries/test_fermi_sign.cpp
  symmetries/test_permutation.cpp
  symmetries/test_permutation_group.cpp
  symmetries/test_generated_group.cpp
  symmetries/test_qn.cpp

  operators/test_opsum.cpp
  operators/test_symmetrize.cpp
  operators/test_non_branching_op.cpp

  blocks/spinhalf/test_spinhalf_matrix.cpp
  blocks/spinhalf/test_spinhalf_apply.cpp

  blocks/spinhalf_symmetric/test_spinhalf_symmetric.cpp
  blocks/spinhalf_symmetric/test_spinhalf_symmetric_matrix.cpp
  blocks/spinhalf_symmetric/test_spinhalf_symmetric_apply.cpp
  blocks/spinhalf_symmetric/test_kitaev_gamma.cpp

  blocks/tj/test_tj_matrix.cpp
  blocks/tj/test_tj_apply.cpp
  blocks/tj/test_tj_raiselower.cpp

  blocks/tj_symmetric/test_tj_symmetric.cpp
  blocks/tj_symmetric/test_tj_symmetric_matrix.cpp
  blocks/tj_symmetric/test_tj_symmetric_apply.cpp

  blocks/electron/test_electron_matrix.cpp
  blocks/electron/test_electron_apply.cpp
  blocks/electron/test_electron_raiselower.cpp

  blocks/electron_symmetric/test_electron_symmetric.cpp
  blocks/electron_symmetric/test_electron_symmetric_matrix.cpp
  blocks/electron_symmetric/test_electron_symmetric_apply.cpp

  algorithms/lanczos/test_eigvals_lanczos.cpp
  algorithms/lanczos/test_eigs_lanczos.cpp
  
  algorithms/lanczos/test_lanczos_pro.cpp
  algorithms/arnoldi/test_arnoldi.cpp
  algorithms/gram_schmidt/test_gram_schmidt.cpp
  algorithms/test_exp_sym_v.cpp
  algorithms/test_norm_estimate.cpp
  algorithms/time_evolution/test_time_evolution.cpp
  algorithms/time_evolution/test_pade.cpp

  states/test_random_state.cpp
  states/test_product_state.cpp
)

set(XDIAG_TESTCASES_SOURCES
  blocks/spinhalf/testcases_spinhalf.cpp
  blocks/tj/testcases_tj.cpp
  blocks/electron/testcases_electron.cpp
)

set(XDIAG_TEST_DISTRIBUTED_SOURCES
  parallel/mpi/test_cdot_distributed.cpp

  basis/spinhalf_distributed/test_basis_sz.cpp
  basis/spinhalf_distributed/test_spinhalf_distributed_basis_iterator.cpp
  basis/tj_distributed/test_basis_np.cpp
  basis/tj_distributed/test_tj_distributed_basis_iterator.cpp

  blocks/spinhalf_distributed/test_spinhalf_distributed.cpp
  blocks/spinhalf_distributed/test_spinhalf_distributed_apply.cpp
  blocks/tj_distributed/test_tj_distributed_apply.cpp

  states/test_product_state_distributed.cpp

  algorithms/time_evolution/test_time_evolution_distributed.cpp
)

add_executable(tests tests.cpp ${XDIAG_TEST_SOURCES} ${XDIAG_TESTCASES_SOURCES})
target_link_libraries(tests PUBLIC ${XDIAG_LIBRARY})
add_test(NAME XdiagTest COMMAND tests)

if(XDIAG_DISTRIBUTED)
  add_executable(tests_distributed tests_distributed.cpp ${XDIAG_TEST_DISTRIBUTED_SOURCES} ${XDIAG_TESTCASES_SOURCES})
  target_link_libraries(tests_distributed PUBLIC ${XDIAG_LIBRARY})
  add_test(NAME XdiagTestDistributed COMMAND tests_distributed)
endif()
