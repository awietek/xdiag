set(XDIAG_SOURCES
  # Utility functions
  utils/format.cpp
  utils/iochecks.cpp
  utils/error.cpp
  utils/logger.cpp
  utils/say_hello.cpp
  utils/read_vectors.cpp
  utils/timing.cpp

  # Input / Output 
  io/read.cpp
  io/file_toml.cpp
  io/toml/file_toml_handler.cpp
  io/toml/value.cpp
  io/toml/std_vector.cpp
  io/toml/arma_vector.cpp
  io/toml/arma_matrix.cpp
  io/toml/operators.cpp

  # Simple math helpers
  math/ipow.cpp
  math/binomial.cpp
  math/scalar.cpp
  math/vector.cpp
  math/matrix.cpp
  math/dot.cpp
  math/norm.cpp
  math/isapprox.cpp

  # Core computational kernels
  kernels/apply.cpp
  kernels/matrix.cpp
  kernels/sparse/coo_matrix.cpp
  kernels/sparse/sparse_build.cpp
  kernels/sparse/csr_matrix.cpp
  kernels/sparse/csc_matrix.cpp
  kernels/sparse/valid.cpp
  kernels/sparse/apply.cpp
  kernels/blocks/spinhalf/kernels.cpp
  kernels/blocks/boson/kernels.cpp
  kernels/blocks/fermion/kernels.cpp
  kernels/blocks/electron/kernels.cpp
  kernels/blocks/tj/kernels.cpp
  kernels/terms/non_branching_op.cpp
  kernels/terms/cdagc_string.cpp

  # Bit(set) functionality
  bits/bitset.cpp
  bits/bitvector.cpp
  bits/bitarray.cpp
  bits/pack_unpack.cpp
  bits/mask_compressor.cpp
  bits/to_string.cpp

  # Combinatorical components
  combinatorics/combinations/enumerate_combinations.cpp
  combinatorics/combinations/combinations.cpp
  combinatorics/combinations/lin_table.cpp
  combinatorics/subsets/subsets.cpp
  combinatorics/bounded_multisets/bounded_multisets.cpp
  combinatorics/bounded_partitions/bounded_partitions.cpp
  combinatorics/bounded_partitions/count_bounded_partitions.cpp
  combinatorics/bounded_partitions/schaefer_table.cpp

  # Bases of the blocks
  basis/basis.cpp
  basis/basis_onthefly.cpp
  basis/basis_electron.cpp
  basis/basis_electron_symmetric.cpp
  basis/basis_tj.cpp
  basis/basis_tj_symmetric.cpp
  basis/basis_sublattice.cpp
  basis/basis_symmetric.cpp

  # Blocks of Hilbert spaces
  blocks/blocks.cpp
  blocks/print_block.cpp
  blocks/spinhalf.cpp
  blocks/boson.cpp
  blocks/fermion.cpp
  blocks/electron.cpp
  blocks/tj.cpp

  # Quantum states
  states/product_state.cpp
  states/random_state.cpp
  states/create_state.cpp
  states/state.cpp
  states/apply.cpp
  states/dot.cpp
  states/inner.cpp
  states/expect.cpp
  states/correlation_matrix.cpp
  states/norm.cpp
  states/fill.cpp

  # utilities for random numbers
  random/random_utils.cpp
  random/hash.cpp

  # Quantum operators
  operators/coeff.cpp
  operators/monomial.cpp
  operators/op.cpp
  operators/opsum.cpp
  operators/types.cpp
  operators/hc.cpp
  operators/valid.cpp
  operators/collect.cpp

  # Permutation symmetries
  symmetries/permutation.cpp
  symmetries/permutation_group.cpp
  symmetries/representation.cpp
  symmetries/representation_set.cpp
  symmetries/cyclic_group.cpp
  symmetries/fermi_sign.cpp
  symmetries/action/site_permutation.cpp
  symmetries/action/isrepresentative.cpp
  symmetries/action/representative.cpp
  symmetries/action/norm.cpp
  symmetries/action/norm_fermionic.cpp
  symmetries/action/site_permutation_sublattice.cpp
  symmetries/action/sublattice_stability.cpp
  symmetries/tables/representative_table.cpp
  symmetries/tables/fermi_table.cpp

  # Algebraic capabilities
  algebra/symmetrize.cpp
  algebra/permute.cpp
  algebra/representation.cpp
  algebra/normal_order.cpp
  algebra/isapprox.cpp
  algebra/ishermitian.cpp
  algebra/algebra.cpp
  algebra/rewrite/rewrite.cpp
  algebra/rewrite/sort_sites.cpp
  algebra/rewrite/spin_rules.cpp
  algebra/rewrite/fermion_rules.cpp
  algebra/rewrite/electron_rules.cpp
  algebra/rewrite/tj_rules.cpp
  algebra/rewrite/matrix_rules.cpp
  algebra/utils/combine_matrix_ops.cpp
  algebra/utils/permute_matrix_op.cpp
  algebra/utils/op_to_matrix_op.cpp
  algebra/utils/swap_pair.cpp
  algebra/utils/replace_pair.cpp
  algebra/utils/check_allowed_types.cpp

  # Numerical linear algebra iterative algorithms
  linalg/lanczos/lanczos_convergence.cpp
  linalg/lanczos/tmatrix.cpp
  linalg/lanczos/eigvals_lanczos.cpp
  linalg/lanczos/eigs_lanczos.cpp
  linalg/lobpcg/eigs_lobpcg.cpp
  linalg/sparse_diag.cpp
  linalg/arnoldi/arnoldi_to_disk.cpp
  linalg/gram_schmidt/gram_schmidt.cpp
  linalg/gram_schmidt/orthogonalize.cpp
  linalg/norm_estimate.cpp
  linalg/time_evolution/time_evolve.cpp
  linalg/time_evolution/imaginary_time_evolve.cpp
  linalg/time_evolution/time_evolve_expokit.cpp
  linalg/time_evolution/evolve_lanczos.cpp
  linalg/time_evolution/expm.cpp
)

set(XDIAG_HDF5_SOURCES
  io/file_h5.cpp
  io/hdf5/file_h5_handler.cpp
  io/hdf5/file_h5_subview.cpp
  io/hdf5/utils.cpp
  io/hdf5/write.cpp
  io/hdf5/types.cpp
)

set(XDIAG_DISTRIBUTED_SOURCES
  # MPI infrastructure
  mpi/allreduce.cpp
  mpi/alltoall.cpp
  mpi/buffer.cpp
  mpi/cdot_distributed.cpp
  mpi/comm_pattern.cpp
  mpi/communicator.cpp
  mpi/datatype.cpp
  mpi/timing_mpi.cpp

  # distributed bases
  basis/distributed/basis_spinhalf_distributed.cpp
  basis/distributed/basis_tj_distributed.cpp
  basis/distributed/basis_electron_distributed.cpp

  # distributed blocks
  blocks/distributed/spinhalf_distributed.cpp
  blocks/distributed/tj_distributed.cpp
  blocks/distributed/electron_distributed.cpp

  # distributed kernels
  kernels/blocks/distributed/apply_distributed.cpp
  kernels/blocks/distributed/spinhalf_distributed/kernels.cpp
  kernels/blocks/distributed/spinhalf_distributed/terms/transpose.cpp
  kernels/blocks/distributed/tj_distributed/kernels.cpp
  kernels/blocks/distributed/electron_distributed/kernels.cpp
)
