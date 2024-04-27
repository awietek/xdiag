set(XDIAG_SOURCES
  utils/iochecks.cpp
  utils/print.cpp
  utils/error.cpp
  utils/logger.cpp
  utils/say_hello.cpp
  utils/read_vectors.cpp
  bits/bitops.cpp
  
  parallel/omp/omp_utils.cpp  
  algebra/algebra.cpp
  algebra/matrix.cpp
  algebra/apply.cpp

  io/args.cpp
  io/args_handler.cpp
  io/file_toml.cpp
  io/toml/file_toml_handler.cpp
  io/toml/toml_conversion.cpp
  io/file_h5.cpp
  io/hdf5/file_h5_handler.cpp
  io/hdf5/utils.cpp
  io/hdf5/write.cpp
  io/hdf5/types.cpp
  
  combinatorics/binomial.cpp
  combinatorics/subsets.cpp
  combinatorics/subsets_index.cpp
  combinatorics/bit_patterns.cpp
  combinatorics/combinations.cpp
  combinatorics/combinations_index.cpp
  combinatorics/combinations_indexing.cpp
  combinatorics/subsets_indexing.cpp
  combinatorics/lin_table.cpp
  combinatorics/fermi_table.cpp
    
  basis/basis.cpp

  basis/spinhalf/basis_sz.cpp
  basis/spinhalf/basis_no_sz.cpp
  basis/spinhalf/basis_symmetric_sz.cpp
  basis/spinhalf/basis_symmetric_no_sz.cpp
  basis/spinhalf/basis_sublattice.cpp

  basis/tj/basis_np.cpp
  basis/tj/basis_symmetric_np.cpp

  basis/electron/basis_np.cpp
  basis/electron/basis_no_np.cpp
  basis/electron/basis_symmetric_np.cpp
  basis/electron/basis_symmetric_no_np.cpp

  blocks/blocks.cpp

  blocks/spinhalf/spinhalf.cpp
  blocks/spinhalf/spinhalf_matrix.cpp
  blocks/spinhalf/spinhalf_apply.cpp
  blocks/spinhalf/compile.cpp
  blocks/spinhalf/qns.cpp

  blocks/electron/electron.cpp
  blocks/electron/electron_matrix.cpp
  blocks/electron/electron_apply.cpp
  blocks/electron/compile.cpp

  blocks/tj/tj.cpp
  blocks/tj/tj_matrix.cpp
  blocks/tj/tj_apply.cpp
  blocks/tj/compile.cpp

  symmetries/qn.cpp
  symmetries/operations/symmetry_operations.cpp	
  symmetries/permutation.cpp
  symmetries/permutation_group.cpp
  symmetries/generated_group.cpp
  symmetries/representation.cpp
  symmetries/operations/fermi_sign.cpp

  symmetries/group_action/group_action.cpp
  symmetries/group_action/group_action_lookup.cpp
  symmetries/group_action/group_action_sublattice.cpp
  symmetries/group_action/sublattice_stability.cpp

  operators/bond.cpp
  operators/bondlist.cpp
  operators/bondlist_handler.cpp
  operators/compiler.cpp
  operators/symmetrized_operator.cpp
  operators/non_branching_bonds.cpp

  states/gpwf.cpp
  states/product_state.cpp
  states/random_state.cpp
  states/state.cpp
  
  random/random_utils.cpp
  random/hash.cpp

  algorithms/lanczos/lanczos_convergence.cpp
  algorithms/lanczos/tmatrix.cpp
  algorithms/lanczos/eigvals_lanczos.cpp
  algorithms/lanczos/eigs_lanczos.cpp
  algorithms/sparse_diag.cpp
  algorithms/arnoldi/arnoldi_to_disk.cpp
  algorithms/gram_schmidt/gram_schmidt.cpp
  algorithms/gram_schmidt/orthogonalize.cpp

  algorithms/norm_estimate.cpp
  algorithms/time_evolution/exp_sym_v.cpp
  algorithms/time_evolution/time_evolution.cpp
  algorithms/time_evolution/pade_matrix_exponential.cpp
)

set(XDIAG_DISTRIBUTED_SOURCES
  parallel/mpi/allreduce.cpp
  parallel/mpi/alltoall.cpp
  parallel/mpi/communicator.cpp
  parallel/mpi/datatype.cpp
  parallel/mpi/cdot_distributed.cpp
  parallel/mpi/timing_mpi.cpp
  parallel/mpi/buffer.cpp

  basis/tj_distributed/basis_np.cpp
  blocks/tj_distributed/tj_distributed.cpp
  blocks/tj_distributed/tj_distributed_apply.cpp
)

set(XDIAG_JULIA_SOURCES
  xdiagjl.cpp
  operators.cpp
  utils.cpp
  blocks.cpp
)
