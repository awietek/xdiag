#pragma once

#if __cplusplus < 201703L
#error "XDiag requires at least C++17"
#endif

#ifdef XDIAG_USE_MPI
#include <mpi.h>
#endif

#ifdef _OPENMP
#include <xdiag/parallel/omp/omp_utils.h>
#include <omp.h>
// #include <xdiag/symmetries/operations/representative_list_omp.h>
#endif

#ifdef XDIAG_USE_HDF5
#define ARMA_USE_HDF5
#endif

#include <xdiag/extern/armadillo/armadillo>
#include <xdiag/extern/clara/clara.hpp>
#include <xdiag/extern/fmt/core.h>

#include <xdiag/common.h>
#include <xdiag/config.h>

#include <xdiag/bits/bitops.h>
#include <xdiag/utils/close.h>
#include <xdiag/utils/iochecks.h>
#include <xdiag/utils/logger.h>
#include <xdiag/utils/precision.h>
#include <xdiag/utils/print.h>
#include <xdiag/utils/print_macro.h>
#include <xdiag/utils/say_hello.h>

#include <xdiag/algebra/algebra.h>
#include <xdiag/algebra/apply.h>
#include <xdiag/algebra/matrix.h>

// Includes for different block types
#include <xdiag/blocks/blocks.h>

#include <xdiag/blocks/spinhalf/compile.h>
#include <xdiag/blocks/spinhalf/spinhalf.h>
#include <xdiag/blocks/spinhalf/spinhalf_apply.h>
#include <xdiag/blocks/spinhalf/spinhalf_matrix.h>

#include <xdiag/blocks/electron/compile.h>
#include <xdiag/blocks/electron/electron.h>
#include <xdiag/blocks/electron/electron_apply.h>
#include <xdiag/blocks/electron/electron_matrix.h>

#include <xdiag/blocks/tj/compile.h>
#include <xdiag/blocks/tj/tj.h>
#include <xdiag/blocks/tj/tj_apply.h>
#include <xdiag/blocks/tj/tj_matrix.h>

#include <xdiag/states/fill.h>
#include <xdiag/states/gpwf.h>
#include <xdiag/states/product_state.h>
#include <xdiag/states/random_state.h>
#include <xdiag/states/state.h>

#include <xdiag/random/hash.h>
#include <xdiag/random/hash_functions.h>
#include <xdiag/random/random_utils.h>

#include <xdiag/combinatorics/binomial.h>
#include <xdiag/combinatorics/bit_patterns.h>
#include <xdiag/combinatorics/combinations.h>
#include <xdiag/combinatorics/combinations_index.h>
#include <xdiag/combinatorics/combinations_indexing.h>
#include <xdiag/combinatorics/fermi_table.h>
#include <xdiag/combinatorics/lin_table.h>
#include <xdiag/combinatorics/subsets.h>
#include <xdiag/combinatorics/subsets_index.h>
#include <xdiag/combinatorics/subsets_indexing.h>

#include <xdiag/basis/basis.h>

#include <xdiag/basis/spinhalf/basis_no_sz.h>
#include <xdiag/basis/spinhalf/basis_sublattice.h>
#include <xdiag/basis/spinhalf/basis_symmetric_no_sz.h>
#include <xdiag/basis/spinhalf/basis_symmetric_sz.h>
#include <xdiag/basis/spinhalf/basis_sz.h>

#include <xdiag/basis/tj/basis_np.h>
#include <xdiag/basis/tj/basis_symmetric_np.h>

#include <xdiag/basis/electron/basis_no_np.h>
#include <xdiag/basis/electron/basis_np.h>
#include <xdiag/basis/electron/basis_symmetric_no_np.h>
#include <xdiag/basis/electron/basis_symmetric_np.h>

#include <xdiag/symmetries/continuous_group.h>
#include <xdiag/symmetries/generated_group.h>
#include <xdiag/symmetries/group_action/group_action.h>
#include <xdiag/symmetries/group_action/group_action_lookup.h>
#include <xdiag/symmetries/group_action/group_action_sublattice.h>
#include <xdiag/symmetries/group_action/sublattice_stability.h>
#include <xdiag/symmetries/operations/fermi_sign.h>
#include <xdiag/symmetries/operations/group_action_operations.h>
#include <xdiag/symmetries/operations/representative_list.h>
#include <xdiag/symmetries/operations/symmetry_operations.h>
#include <xdiag/symmetries/permutation.h>
#include <xdiag/symmetries/permutation_group.h>
#include <xdiag/symmetries/qn.h>
#include <xdiag/symmetries/representation.h>

#include <xdiag/operators/bond.h>
#include <xdiag/operators/bondlist.h>
#include <xdiag/operators/bondlist_handler.h>
#include <xdiag/operators/compiler.h>
#include <xdiag/operators/non_branching_bonds.h>
#include <xdiag/operators/symmetrized_operator.h>


#include <xdiag/algorithms/norm_estimate.h>
#include <xdiag/algorithms/sparse_diag.h>

#include <xdiag/algorithms/lanczos/lanczos.h>
#include <xdiag/algorithms/lanczos/eigs_lanczos.h>
#include <xdiag/algorithms/lanczos/eigvals_lanczos.h>
#include <xdiag/algorithms/lanczos/lanczos_convergence.h>
#include <xdiag/algorithms/lanczos/tmatrix.h>
#include <xdiag/algorithms/lanczos/lanczos_pro.h>

#include <xdiag/algorithms/arnoldi/arnoldi.h>
#include <xdiag/algorithms/arnoldi/arnoldi_to_disk.h>
#include <xdiag/algorithms/gram_schmidt/gram_schmidt.h>
#include <xdiag/algorithms/gram_schmidt/orthogonalize.h>

#include <xdiag/algorithms/time_evolution/exp_sym_v.h>
#include <xdiag/algorithms/time_evolution/pade_matrix_exponential.h>
#include <xdiag/algorithms/time_evolution/time_evolution.h>
#include <xdiag/algorithms/time_evolution/zahexpv.h>

#include <xdiag/io/args.h>
#include <xdiag/io/file_h5.h>
#include <xdiag/io/file_toml.h>

#ifdef XDIAG_ENABLE_MPI
#include <xdiag/mpi/allreduce.h>
#include <xdiag/mpi/alltoall.h>
#include <xdiag/mpi/buffer.h>
#include <xdiag/mpi/datatype.h>
#include <xdiag/mpi/dot_mpi.h>
#include <xdiag/mpi/logger_mpi.h>
#include <xdiag/mpi/timing_mpi.h>
#include <xdiag/utils/print_mpi.h>

#include <xdiag/basis/spinhalf_mpi/spinhalf_mpi_basis_sz.h>
#include <xdiag/blocks/utils/block_utils_mpi.h>

#include <xdiag/blocks/spinhalf_mpi/spinhalf_mpi.h>
#include <xdiag/blocks/spinhalf_mpi/spinhalf_mpi_apply.h>
#include <xdiag/blocks/spinhalf_mpi/spinhalf_mpi_fill.h>
#include <xdiag/blocks/spinhalf_mpi/terms/get_prefix_postfix_mixed_bonds.h>

#include <xdiag/blocks/electron_mpi/electron_mpi.h>
#include <xdiag/blocks/electron_mpi/electron_mpi_apply.h>
#endif
