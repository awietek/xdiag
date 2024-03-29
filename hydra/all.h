#pragma once

#if __cplusplus < 201703L
#error "Hydra requires at least C++17"
#endif

#ifdef HYDRA_USE_MPI
#include <mpi.h>
#endif

#ifdef _OPENMP
#include <hydra/parallel/omp/omp_utils.h>
#include <omp.h>
// #include <hydra/symmetries/operations/representative_list_omp.h>
#endif

#ifdef HYDRA_USE_HDF5
#define ARMA_USE_HDF5
#endif

#include <hydra/extern/armadillo/armadillo>
#include <hydra/extern/clara/clara.hpp>
#include <hydra/extern/fmt/core.h>

#include <hydra/common.h>
#include <hydra/config.h>

#include <hydra/bits/bitops.h>
#include <hydra/utils/close.h>
#include <hydra/utils/iochecks.h>
#include <hydra/utils/logger.h>
#include <hydra/utils/precision.h>
#include <hydra/utils/print.h>
#include <hydra/utils/print_macro.h>
#include <hydra/utils/say_hello.h>

#include <hydra/algebra/algebra.h>
#include <hydra/algebra/apply.h>
#include <hydra/algebra/matrix.h>

// Includes for different block types
#include <hydra/blocks/blocks.h>

#include <hydra/blocks/spinhalf/compile.h>
#include <hydra/blocks/spinhalf/spinhalf.h>
#include <hydra/blocks/spinhalf/spinhalf_apply.h>
#include <hydra/blocks/spinhalf/spinhalf_matrix.h>

#include <hydra/blocks/electron/compile.h>
#include <hydra/blocks/electron/electron.h>
#include <hydra/blocks/electron/electron_apply.h>
#include <hydra/blocks/electron/electron_matrix.h>

#include <hydra/blocks/tj/compile.h>
#include <hydra/blocks/tj/tj.h>
#include <hydra/blocks/tj/tj_apply.h>
#include <hydra/blocks/tj/tj_matrix.h>

#include <hydra/states/fill.h>
#include <hydra/states/gpwf.h>
#include <hydra/states/product_state.h>
#include <hydra/states/random_state.h>
#include <hydra/states/state.h>

#include <hydra/random/hash.h>
#include <hydra/random/hash_functions.h>
#include <hydra/random/random_utils.h>

#include <hydra/combinatorics/binomial.h>
#include <hydra/combinatorics/bit_patterns.h>
#include <hydra/combinatorics/combinations.h>
#include <hydra/combinatorics/combinations_index.h>
#include <hydra/combinatorics/combinations_indexing.h>
#include <hydra/combinatorics/fermi_table.h>
#include <hydra/combinatorics/lin_table.h>
#include <hydra/combinatorics/subsets.h>
#include <hydra/combinatorics/subsets_index.h>
#include <hydra/combinatorics/subsets_indexing.h>

#include <hydra/basis/basis.h>

#include <hydra/basis/spinhalf/basis_no_sz.h>
#include <hydra/basis/spinhalf/basis_sublattice.h>
#include <hydra/basis/spinhalf/basis_symmetric_no_sz.h>
#include <hydra/basis/spinhalf/basis_symmetric_sz.h>
#include <hydra/basis/spinhalf/basis_sz.h>

#include <hydra/basis/tj/basis_np.h>
#include <hydra/basis/tj/basis_symmetric_np.h>

#include <hydra/basis/electron/basis_no_np.h>
#include <hydra/basis/electron/basis_np.h>
#include <hydra/basis/electron/basis_symmetric_no_np.h>
#include <hydra/basis/electron/basis_symmetric_np.h>

#include <hydra/symmetries/continuous_group.h>
#include <hydra/symmetries/generated_group.h>
#include <hydra/symmetries/group_action/group_action.h>
#include <hydra/symmetries/group_action/group_action_lookup.h>
#include <hydra/symmetries/group_action/group_action_sublattice.h>
#include <hydra/symmetries/group_action/sublattice_stability.h>
#include <hydra/symmetries/operations/fermi_sign.h>
#include <hydra/symmetries/operations/group_action_operations.h>
#include <hydra/symmetries/operations/representative_list.h>
#include <hydra/symmetries/operations/symmetry_operations.h>
#include <hydra/symmetries/permutation.h>
#include <hydra/symmetries/permutation_group.h>
#include <hydra/symmetries/qn.h>
#include <hydra/symmetries/representation.h>

#include <hydra/operators/bond.h>
#include <hydra/operators/bondlist.h>
#include <hydra/operators/bondlist_handler.h>
#include <hydra/operators/compiler.h>
#include <hydra/operators/non_branching_bonds.h>
#include <hydra/operators/symmetrized_operator.h>


#include <hydra/algorithms/norm_estimate.h>
#include <hydra/algorithms/sparse_diag.h>

#include <hydra/algorithms/lanczos/lanczos.h>
#include <hydra/algorithms/lanczos/eigs_lanczos.h>
#include <hydra/algorithms/lanczos/eigvals_lanczos.h>
#include <hydra/algorithms/lanczos/lanczos_convergence.h>
#include <hydra/algorithms/lanczos/tmatrix.h>
#include <hydra/algorithms/lanczos/lanczos_pro.h>

#include <hydra/algorithms/arnoldi/arnoldi.h>
#include <hydra/algorithms/arnoldi/arnoldi_to_disk.h>
#include <hydra/algorithms/gram_schmidt/gram_schmidt.h>
#include <hydra/algorithms/gram_schmidt/orthogonalize.h>

#include <hydra/algorithms/time_evolution/exp_sym_v.h>
#include <hydra/algorithms/time_evolution/pade_matrix_exponential.h>
#include <hydra/algorithms/time_evolution/time_evolution.h>
#include <hydra/algorithms/time_evolution/zahexpv.h>

#include <hydra/io/args.h>
#include <hydra/io/file_h5.h>
#include <hydra/io/file_toml.h>

#ifdef HYDRA_ENABLE_MPI
#include <hydra/mpi/allreduce.h>
#include <hydra/mpi/alltoall.h>
#include <hydra/mpi/buffer.h>
#include <hydra/mpi/datatype.h>
#include <hydra/mpi/dot_mpi.h>
#include <hydra/mpi/logger_mpi.h>
#include <hydra/mpi/timing_mpi.h>
#include <hydra/utils/print_mpi.h>

#include <hydra/basis/spinhalf_mpi/spinhalf_mpi_basis_sz.h>
#include <hydra/blocks/utils/block_utils_mpi.h>

#include <hydra/blocks/spinhalf_mpi/spinhalf_mpi.h>
#include <hydra/blocks/spinhalf_mpi/spinhalf_mpi_apply.h>
#include <hydra/blocks/spinhalf_mpi/spinhalf_mpi_fill.h>
#include <hydra/blocks/spinhalf_mpi/terms/get_prefix_postfix_mixed_bonds.h>

#include <hydra/blocks/electron_mpi/electron_mpi.h>
#include <hydra/blocks/electron_mpi/electron_mpi_apply.h>
#endif
