set(HYDRA_SOURCES
  hydra/utils/iochecks.cpp
  hydra/utils/print.cpp
  hydra/bitops/bitops.cpp
  
  hydra/parallel/omp/omp_utils.cpp
  
  hydra/algebra/algebra.cpp
  hydra/algebra/matrix.cpp
  
  hydra/io/file_toml.cpp
  hydra/io/toml/file_toml_handler.cpp
  hydra/io/toml/toml_conversion.cpp
  hydra/io/file_h5.cpp
  hydra/io/hdf5/file_h5_handler.cpp
  hydra/io/hdf5/utils.cpp
  hydra/io/hdf5/write.cpp
  hydra/io/hdf5/types.cpp
  
  hydra/combinatorics/binomial.cpp
  hydra/combinatorics/subsets.cpp
  hydra/combinatorics/subsets_index.cpp
  hydra/combinatorics/bit_patterns.cpp
  hydra/combinatorics/combinations.cpp
  hydra/combinatorics/combinations_index.cpp
  
  hydra/indexing/indexing_variants.cpp
  hydra/indexing/lin_table.cpp
  hydra/indexing/fermi_table.cpp
  hydra/indexing/combinations_indexing.cpp
  hydra/indexing/subsets_indexing.cpp
  
  hydra/indexing/spinhalf/indexing_sz.cpp
  hydra/indexing/spinhalf/indexing_no_sz.cpp
  hydra/indexing/spinhalf/indexing_symmetric_sz.cpp
  hydra/indexing/spinhalf/indexing_symmetric_no_sz.cpp
  hydra/indexing/spinhalf/indexing_sublattice.cpp
  hydra/indexing/spinhalf/symmetric_iterator.cpp

  hydra/indexing/tj/indexing_np.cpp
  hydra/indexing/tj/indexing_symmetric_np.cpp

  hydra/indexing/electron/indexing_np.cpp
  hydra/indexing/electron/indexing_no_np.cpp
  hydra/indexing/electron/indexing_symmetric_np.cpp
  hydra/indexing/electron/indexing_symmetric_no_np.cpp

  hydra/blocks/blocks.cpp
  hydra/blocks/utils/block_utils.cpp
  hydra/blocks/spinhalf/spinhalf.cpp
  hydra/blocks/spinhalf/spinhalf_matrix.cpp
  hydra/blocks/spinhalf/spinhalf_apply.cpp
  hydra/blocks/spinhalf/terms/compile.cpp
  hydra/blocks/spinhalf/terms/qns.cpp

  hydra/blocks/electron/electron.cpp
  hydra/blocks/electron/electron_matrix.cpp
  hydra/blocks/electron/electron_apply.cpp
  hydra/blocks/electron/terms/compile.cpp

  hydra/blocks/tj/tj.cpp
  hydra/blocks/tj/tj_matrix.cpp
  hydra/blocks/tj/tj_apply.cpp
  hydra/blocks/tj/terms/compile.cpp

  hydra/symmetries/operations/symmetry_operations.cpp
  hydra/symmetries/permutation.cpp
  hydra/symmetries/permutation_group.cpp
  hydra/symmetries/generated_group.cpp
  hydra/symmetries/representation.cpp
  hydra/symmetries/operations/fermi_sign.cpp

  hydra/symmetries/group_action/group_action.cpp
  hydra/symmetries/group_action/group_action_lookup.cpp
  hydra/symmetries/group_action/group_action_sublattice.cpp
  hydra/symmetries/group_action/sublattice_stability.cpp

  hydra/operators/bond.cpp
  hydra/operators/bondlist.cpp
  hydra/operators/bondlist_handler.cpp
  hydra/operators/compiler.cpp
  hydra/operators/symmetrized_operator.cpp
  hydra/operators/non_branching_bonds.cpp

  hydra/states/gpwf_spinhalf.cpp
  hydra/states/product_state.cpp
  hydra/states/random_state.cpp
  hydra/states/state.cpp

  hydra/random/random_utils.cpp
  hydra/random/hashes.cpp

  hydra/algorithms/lanczos/lanczos_convergence.cpp
  hydra/algorithms/lanczos/tmatrix.cpp
  hydra/algorithms/lanczos/lanczos_eigenvector.cpp
  hydra/algorithms/lanczos/lanczos_eigenvalues.cpp
  hydra/algorithms/sparse_diag.cpp
  hydra/algorithms/norm_estimate.cpp
  hydra/algorithms/time_evolution/zahexpv.cpp
  hydra/algorithms/time_evolution/time_evolution.cpp
  hydra/algorithms/time_evolution/pade_matrix_exponential.cpp
)
