sources+= hydra/utils/iochecks.cpp
sources+= hydra/utils/print.cpp
sources+= hydra/bitops/bitops.cpp

sources+= hydra/parallel/omp/omp_utils.cpp

sources+= hydra/combinatorics/binomial.cpp
sources+= hydra/combinatorics/subsets.cpp
sources+= hydra/combinatorics/subsets_index.cpp
sources+= hydra/combinatorics/bit_patterns.cpp
sources+= hydra/combinatorics/combinations.cpp
sources+= hydra/combinatorics/combinations_index.cpp

sources+= hydra/indexing/lin_table.cpp
sources+= hydra/indexing/fermi_table.cpp
sources+= hydra/indexing/combinations_indexing.cpp
sources+= hydra/indexing/subsets_indexing.cpp

sources+= hydra/indexing/spinhalf/indexing_sz.cpp
sources+= hydra/indexing/spinhalf/indexing_no_sz.cpp
sources+= hydra/indexing/spinhalf/indexing_symmetric_sz.cpp
sources+= hydra/indexing/spinhalf/indexing_symmetric_no_sz.cpp
sources+= hydra/indexing/spinhalf/indexing_sublattice.cpp
sources+= hydra/indexing/spinhalf/symmetric_iterator.cpp

sources+= hydra/indexing/tj/indexing_np.cpp
sources+= hydra/indexing/tj/indexing_symmetric_np.cpp

sources+= hydra/indexing/electron/indexing_np.cpp
sources+= hydra/indexing/electron/indexing_no_np.cpp
sources+= hydra/indexing/electron/indexing_symmetric_np.cpp
sources+= hydra/indexing/electron/indexing_symmetric_no_np.cpp

sources+= hydra/blocks/utils/block_utils.cpp

sources+= hydra/blocks/spinhalf/spinhalf.cpp
sources+= hydra/blocks/spinhalf/spinhalf_matrix.cpp
sources+= hydra/blocks/spinhalf/spinhalf_apply.cpp
sources+= hydra/blocks/spinhalf/terms/compile.cpp
sources+= hydra/blocks/spinhalf/terms/qns.cpp

sources+= hydra/blocks/electron/electron.cpp
sources+= hydra/blocks/electron/electron_matrix.cpp
sources+= hydra/blocks/electron/electron_apply.cpp
sources+= hydra/blocks/electron/terms/compile.cpp

sources+= hydra/blocks/tj/tj.cpp
sources+= hydra/blocks/tj/tj_matrix.cpp
sources+= hydra/blocks/tj/tj_apply.cpp
sources+= hydra/blocks/tj/terms/compile.cpp

sources+= hydra/symmetries/operations/symmetry_operations.cpp
sources+= hydra/symmetries/permutation.cpp
sources+= hydra/symmetries/permutation_group.cpp
sources+= hydra/symmetries/representation.cpp
sources+= hydra/symmetries/operations/fermi_sign.cpp

sources+= hydra/symmetries/group_action/group_action.cpp
sources+= hydra/symmetries/group_action/group_action_lookup.cpp
sources+= hydra/symmetries/group_action/group_action_sublattice.cpp
sources+= hydra/symmetries/group_action/sublattice_stability.cpp

sources+= hydra/operators/bond.cpp
sources+= hydra/operators/bondlist.cpp
sources+= hydra/operators/compiler.cpp
sources+= hydra/operators/symmetric_operator.cpp
sources+= hydra/operators/non_branching_bonds.cpp

sources+= hydra/states/gpwf_spinhalf.cpp

sources+= hydra/random/random_utils.cpp
sources+= hydra/random/hashes.cpp

sources+=hydra/algorithms/lanczos/lanczos_convergence.cpp
sources+=hydra/algorithms/lanczos/tmatrix.cpp

mpisources+= hydra/utils/print_mpi.cpp

mpisources+= hydra/mpi/datatype.cpp
mpisources+= hydra/mpi/allreduce.cpp
mpisources+= hydra/mpi/alltoall.cpp
mpisources+= hydra/mpi/timing_mpi.cpp
mpisources+= hydra/mpi/dot_mpi.cpp
mpisources+= hydra/mpi/communicator.cpp
mpisources+= hydra/mpi/buffer.cpp

mpisources+= hydra/indexing/spinhalf_mpi/spinhalf_mpi_indexing_sz.cpp

mpisources+= hydra/blocks/utils/block_utils_mpi.cpp
mpisources+= hydra/blocks/spinhalf_mpi/terms/get_prefix_postfix_mixed_bonds.cpp
mpisources+= hydra/blocks/spinhalf_mpi/spinhalf_mpi.cpp
mpisources+= hydra/blocks/spinhalf_mpi/spinhalf_mpi_apply.cpp
mpisources+= hydra/blocks/electron_mpi/electron_mpi.cpp
mpisources+= hydra/blocks/electron_mpi/electron_mpi_apply.cpp


testsources+= test/tests.cpp
testsources+= test/bitops/test_bitops.cpp

testsources+= test/combinatorics/test_binomial.cpp
testsources+= test/combinatorics/test_subsets.cpp
testsources+= test/combinatorics/test_bit_patterns.cpp
testsources+= test/combinatorics/test_combinations.cpp
testsources+= test/combinatorics/test_combinations_index.cpp

testsources+= test/indexing/test_lin_table.cpp
testsources+= test/indexing/test_fermi_table.cpp
testsources+= test/indexing/spinhalf/test_spinhalf_indexing_sublattice.cpp
testsources+= test/indexing/spinhalf/test_spinhalf_indexing.cpp

testsources+= test/symmetries/group_action/test_group_action.cpp
testsources+= test/symmetries/group_action/test_group_action_lookup.cpp
testsources+= test/symmetries/group_action/test_group_action_sublattice.cpp
testsources+= test/symmetries/test_fermi_sign.cpp
testsources+= test/symmetries/test_symmetry_operations.cpp
testsources+= test/symmetries/test_permutation.cpp
testsources+= test/symmetries/test_permutation_group.cpp

testsources+= test/operators/test_symmetric_operator.cpp
testsources+= test/operators/test_non_branching_bonds.cpp

testsources+= test/blocks/spinhalf/testcases_spinhalf.cpp
testsources+= test/blocks/spinhalf/test_spinhalf_matrix.cpp
testsources+= test/blocks/spinhalf/test_spinhalf_apply.cpp

testsources+= test/blocks/spinhalf_symmetric/test_spinhalf_symmetric.cpp
testsources+= test/blocks/spinhalf_symmetric/test_spinhalf_symmetric_matrix.cpp
testsources+= test/blocks/spinhalf_symmetric/test_spinhalf_symmetric_apply.cpp

testsources+= test/blocks/tj/testcases_tj.cpp
testsources+= test/blocks/tj/test_tj_utils.cpp
testsources+= test/blocks/tj/test_tj_matrix.cpp
testsources+= test/blocks/tj/test_tj_apply.cpp

testsources+= test/blocks/tj_symmetric/test_tj_symmetric.cpp
testsources+= test/blocks/tj_symmetric/test_tj_symmetric_matrix.cpp
testsources+= test/blocks/tj_symmetric/test_tj_symmetric_apply.cpp

testsources+= test/blocks/electron/testcases_electron.cpp
testsources+= test/blocks/electron/test_electron_matrix.cpp
testsources+= test/blocks/electron/test_electron_apply.cpp

testsources+= test/blocks/electron_symmetric/test_electron_symmetric.cpp
testsources+= test/blocks/electron_symmetric/test_electron_symmetric_matrix.cpp
testsources+= test/blocks/electron_symmetric/test_electron_symmetric_apply.cpp

testsources+= test/algorithms/lanczos/test_lanczos_eigenvalues.cpp
testsources+= test/algorithms/lanczos/test_lanczos_eigenvector.cpp

testsources+= test/states/test_random_state.cpp

# testmpisources+= testmpi/tests.cpp
# testmpisources+= testmpi/mpi/test_dot_mpi.cpp
# testmpisources+= test/blocks/spinhalf/testcases_spinhalf.cpp
# testmpisources+= testmpi/blocks/spinhalf_mpi/test_spinhalf_mpi.cpp
# testmpisources+= testmpi/blocks/spinhalf_mpi/test_spinhalf_mpi_apply.cpp
