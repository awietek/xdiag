set(XDIAG_SOURCES
  utils/iochecks.cpp
  utils/error.cpp
  utils/logger.cpp
  utils/say_hello.cpp
  utils/read_vectors.cpp
  utils/type_string.cpp

  utils/scalar.cpp
  utils/vector.cpp
  utils/matrix.cpp

  bits/bitops.cpp
  
  parallel/omp/omp_utils.cpp  
  algebra/algebra.cpp
  algebra/matrix.cpp
  algebra/apply.cpp

  io/read.cpp
  io/file_toml.cpp
  io/file_h5.cpp
  io/toml/file_toml_handler.cpp
  # io/toml/toml_conversion.cpp
  io/toml/value.cpp
  io/toml/std_vector.cpp
  io/toml/arma_vector.cpp
  io/toml/arma_matrix.cpp
  io/toml/operators.cpp
  io/hdf5/file_h5_handler.cpp
  io/hdf5/file_h5_subview.cpp
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

  basis/spinhalf/basis_spinhalf.cpp
  basis/spinhalf/basis_sz.cpp
  basis/spinhalf/basis_no_sz.cpp
  basis/spinhalf/basis_symmetric_sz.cpp
  basis/spinhalf/basis_symmetric_no_sz.cpp
  basis/spinhalf/basis_sublattice.cpp
  basis/spinhalf/apply/dispatch_matrix.cpp
  basis/spinhalf/apply/dispatch_apply.cpp

  basis/tj/basis_tj.cpp
  basis/tj/basis_np.cpp
  basis/tj/basis_symmetric_np.cpp
  basis/tj/apply/dispatch_matrix.cpp
  basis/tj/apply/dispatch_apply.cpp
  
  basis/electron/basis_electron.cpp
  basis/electron/basis_np.cpp
  basis/electron/basis_no_np.cpp
  basis/electron/basis_symmetric_np.cpp
  basis/electron/basis_symmetric_no_np.cpp
  basis/electron/apply/dispatch_matrix.cpp
  basis/electron/apply/dispatch_apply.cpp
  
  blocks/blocks.cpp
  blocks/spinhalf.cpp
  blocks/electron.cpp
  blocks/tj.cpp

  symmetries/permutation.cpp
  symmetries/permutation_group.cpp
  symmetries/representation.cpp
  symmetries/operations/symmetry_operations.cpp	
  symmetries/operations/fermi_sign.cpp
  symmetries/group_action/group_action.cpp
  symmetries/group_action/group_action_lookup.cpp
  symmetries/group_action/group_action_sublattice.cpp
  symmetries/group_action/sublattice_stability.cpp

  operators/coupling.cpp
  operators/op.cpp
  operators/opsum.cpp
  operators/logic/compilation.cpp
  operators/logic/valid.cpp
  operators/logic/types.cpp
  operators/logic/symmetrize.cpp
  operators/logic/real.cpp
  operators/logic/hc.cpp
  operators/logic/isapprox.cpp
  operators/logic/permute.cpp
  operators/logic/qns.cpp
  operators/logic/non_branching_op.cpp
  operators/logic/order.cpp
  
  states/gpwf.cpp
  states/product_state.cpp
  states/random_state.cpp
  states/state.cpp
  states/fill.cpp
  states/create_state.cpp
  
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
  parallel/mpi/comm_pattern.cpp
  parallel/mpi/datatype.cpp
  parallel/mpi/cdot_distributed.cpp
  parallel/mpi/timing_mpi.cpp
  parallel/mpi/buffer.cpp

  basis/spinhalf_distributed/basis_spinhalf_distributed.cpp
  basis/spinhalf_distributed/basis_sz.cpp
  basis/spinhalf_distributed/transpose.cpp
  basis/spinhalf_distributed/apply/dispatch_apply.cpp
  basis/spinhalf_distributed/apply/apply_terms.cpp
  
  basis/tj_distributed/basis_tj_distributed.cpp
  basis/tj_distributed/basis_np.cpp
  basis/tj_distributed/apply/dispatch_apply.cpp

  blocks/spinhalf_distributed.cpp
  blocks/tj_distributed.cpp
)

set(XDIAG_JULIA_SOURCES
  xdiagjl.cpp
  utils/armadillo.cpp
  algebra/matrix.cpp
  algebra/apply.cpp
  algebra/algebra.cpp
  algorithms/sparse_diag.cpp
  algorithms/lanczos/eigs_lanczos.cpp
  algorithms/lanczos/eigvals_lanczos.cpp
  blocks/spinhalf.cpp
  blocks/tj.cpp
  blocks/electron.cpp
  operators/coupling.cpp
  operators/op.cpp
  operators/opsum.cpp
  operators/symmetrize.cpp
  states/create_state.cpp
  states/fill.cpp
  states/random_state.cpp
  states/product_state.cpp
  states/gpwf.cpp
  states/state.cpp
  utils/utils.cpp
  symmetries/permutation.cpp
  symmetries/permutation_group.cpp
  symmetries/representation.cpp
)
