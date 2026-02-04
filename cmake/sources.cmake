set(XDIAG_SOURCES
  utils/iochecks.cpp
  utils/error.cpp
  utils/logger.cpp
  utils/say_hello.cpp
  utils/read_vectors.cpp
  utils/type_string.cpp
  utils/split.cpp
  utils/arma_to_cx.cpp
  utils/scalar.cpp
  utils/vector.cpp
  utils/matrix.cpp

  bits/bitops.cpp
  bits/bitset.cpp
  bits/bitvector.cpp

  # combinatorics/binomial.cpp
  # combinatorics/bit_patterns.cpp
  # combinatorics/combinations/combinations.cpp
  # combinatorics/combinations/combinations_table.cpp
  # combinatorics/combinations/lin_table.cpp
  # combinatorics/subsets/subsets.cpp
  # combinatorics/subsets/subsets_table.cpp
    
  # parallel/omp/omp_utils.cpp  
  # algebra/algebra.cpp
  # algebra/matrix.cpp
  # algebra/apply.cpp
  # algebra/apply_dispatch.cpp
  # algebra/isapprox.cpp
  
  # algebra/sparse/apply.cpp
  
  # algebra/sparse/coo_matrix.cpp
  # algebra/sparse/coo_matrix_generate.cpp
  # algebra/sparse/instantiations/coo_matrix_nnz_double.cpp
  # algebra/sparse/instantiations/coo_matrix_nnz_complex.cpp
  # algebra/sparse/instantiations/coo_matrix_fill_int32_double.cpp
  # algebra/sparse/instantiations/coo_matrix_fill_int32_complex.cpp
  # algebra/sparse/instantiations/coo_matrix_fill_int64_double.cpp
  # algebra/sparse/instantiations/coo_matrix_fill_int64_complex.cpp
  
  # algebra/sparse/csr_matrix.cpp
  # algebra/sparse/csr_matrix_generate.cpp
  # algebra/sparse/instantiations/csr_matrix_nnz_int32_double.cpp
  # algebra/sparse/instantiations/csr_matrix_nnz_int32_complex.cpp
  # algebra/sparse/instantiations/csr_matrix_nnz_int64_double.cpp
  # algebra/sparse/instantiations/csr_matrix_nnz_int64_complex.cpp
  # algebra/sparse/instantiations/csr_matrix_fill_int32_double.cpp
  # algebra/sparse/instantiations/csr_matrix_fill_int32_complex.cpp
  # algebra/sparse/instantiations/csr_matrix_fill_int64_double.cpp
  # algebra/sparse/instantiations/csr_matrix_fill_int64_complex.cpp
  
  # algebra/sparse/csc_matrix.cpp
  
  # io/read.cpp
  # io/file_toml.cpp
  # io/file_h5.cpp
  # io/toml/file_toml_handler.cpp
  # io/toml/value.cpp
  # io/toml/std_vector.cpp
  # io/toml/arma_vector.cpp
  # io/toml/arma_matrix.cpp
  # io/toml/operators.cpp
  # io/hdf5/file_h5_handler.cpp
  # io/hdf5/file_h5_subview.cpp
  # io/hdf5/utils.cpp
  # io/hdf5/write.cpp
  # io/hdf5/types.cpp
  
  # basis/spinhalf/basis_spinhalf.cpp
  # basis/spinhalf/basis_sz.cpp
  # basis/spinhalf/basis_no_sz.cpp
  # basis/spinhalf/basis_symmetric_sz.cpp
  # basis/spinhalf/basis_symmetric_no_sz.cpp
  # basis/spinhalf/basis_sublattice.cpp

  # basis/tj/basis_tj.cpp
  # basis/tj/basis_np.cpp
  # basis/tj/basis_symmetric_np.cpp

  # basis/electron/basis_electron.cpp
  # basis/electron/basis_np.cpp
  # basis/electron/basis_no_np.cpp
  # basis/electron/basis_symmetric_np.cpp
  # basis/electron/basis_symmetric_no_np.cpp
  
  blocks/blocks.cpp
  blocks/spinhalf.cpp
  blocks/electron.cpp
  blocks/tj.cpp

  # symmetries/permutation.cpp
  # symmetries/permutation_group.cpp
  # symmetries/representation.cpp
  # symmetries/operations/symmetry_operations.cpp	
  # symmetries/operations/representative_lookup_table.cpp	
  # symmetries/operations/fermi_sign.cpp
  # symmetries/group_action/group_action.cpp
  # symmetries/group_action/group_action_lookup.cpp
  # symmetries/group_action/group_action_sublattice.cpp
  # symmetries/group_action/sublattice_stability.cpp

  # operators/coupling.cpp
  # operators/op.cpp
  # operators/opsum.cpp
  # operators/logic/compilation.cpp
  # operators/logic/valid.cpp
  # operators/logic/types.cpp
  # operators/logic/symmetrize.cpp
  # operators/logic/real.cpp
  # operators/logic/hc.cpp
  # operators/logic/isapprox.cpp
  # operators/logic/permute.cpp
  # operators/logic/qns.cpp
  # operators/logic/non_branching_op.cpp
  # operators/logic/order.cpp
  # operators/logic/block.cpp
  
  # states/gpwf.cpp
  # states/product_state.cpp
  # states/random_state.cpp
  # states/state.cpp
  # states/fill.cpp
  # states/create_state.cpp
  
  # random/random_utils.cpp
  # random/hash.cpp

  # algorithms/lanczos/lanczos_convergence.cpp
  # algorithms/lanczos/tmatrix.cpp
  # algorithms/lanczos/eigvals_lanczos.cpp
  # algorithms/lanczos/eigs_lanczos.cpp
  # algorithms/sparse_diag.cpp
  # algorithms/arnoldi/arnoldi_to_disk.cpp
  # algorithms/gram_schmidt/gram_schmidt.cpp
  # algorithms/gram_schmidt/orthogonalize.cpp

  # algorithms/norm_estimate.cpp
  # algorithms/time_evolution/time_evolve.cpp
  # algorithms/time_evolution/imaginary_time_evolve.cpp
  # algorithms/time_evolution/time_evolve_expokit.cpp
  # algorithms/time_evolution/evolve_lanczos.cpp
  # algorithms/time_evolution/expm.cpp
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
  basis/spinhalf_distributed/apply/apply_terms.cpp
  
  basis/tj_distributed/basis_tj_distributed.cpp
  basis/tj_distributed/basis_np.cpp

  basis/electron_distributed/basis_electron_distributed.cpp
  basis/electron_distributed/basis_np.cpp
  
  blocks/spinhalf_distributed.cpp
  blocks/tj_distributed.cpp
  blocks/electron_distributed.cpp
)

set(XDIAG_JULIA_SOURCES
  xdiagjl.cpp

  algebra/matrix.cpp
  algebra/apply.cpp
  algebra/algebra.cpp
  algebra/sparse/coo_matrix.cpp
  algebra/sparse/csr_matrix.cpp
  algebra/sparse/csc_matrix.cpp
  algebra/sparse/apply.cpp

  algorithms/sparse_diag.cpp
  algorithms/lanczos/eigs_lanczos.cpp
  algorithms/lanczos/eigvals_lanczos.cpp
  algorithms/time_evolution/time_evolve.cpp
  algorithms/time_evolution/imaginary_time_evolve.cpp
  algorithms/time_evolution/time_evolve_expokit.cpp  
  algorithms/time_evolution/evolve_lanczos.cpp

  blocks/spinhalf.cpp
  blocks/tj.cpp
  blocks/electron.cpp

  io/file_toml.cpp
  io/read.cpp
  
  operators/op.cpp
  operators/opsum.cpp
  operators/symmetrize.cpp
  operators/hc.cpp
  operators/block.cpp
  
  states/create_state.cpp
  states/fill.cpp
  states/random_state.cpp
  states/product_state.cpp
  states/gpwf.cpp
  states/state.cpp

  utils/utils.cpp
  utils/armadillo.cpp

  symmetries/permutation.cpp
  symmetries/permutation_group.cpp
  symmetries/representation.cpp
)
