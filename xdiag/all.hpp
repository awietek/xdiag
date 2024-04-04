#pragma once

#if __cplusplus < 201703L
#error "XDiag requires at least C++17"
#endif

#ifdef XDIAG_USE_MPI
#include <mpi.h>
#endif

#ifdef _OPENMP
#include <xdiag/parallel/omp/omp_utils.hpp>
#include <omp.h>
// #include <xdiag/symmetries/operations/representative_list_omp.hpp>
#endif

#ifdef XDIAG_USE_HDF5
#define ARMA_USE_HDF5
#endif

#include <xdiag/extern/armadillo/armadillo>
#include <xdiag/extern/clara/clara.hpp>
#include <xdiag/extern/fmt/core.hpp>

#include <xdiag/common.hpp>
#include <xdiag/config.hpp>

#include <xdiag/bits/bitops.hpp>
#include <xdiag/utils/close.hpp>
#include <xdiag/utils/iochecks.hpp>
#include <xdiag/utils/logger.hpp>
#include <xdiag/utils/precision.hpp>
#include <xdiag/utils/print.hpp>
#include <xdiag/utils/print_macro.hpp>
#include <xdiag/utils/say_hello.hpp>

#include <xdiag/algebra/algebra.hpp>
#include <xdiag/algebra/apply.hpp>
#include <xdiag/algebra/matrix.hpp>

// Includes for different block types
#include <xdiag/blocks/blocks.hpp>

#include <xdiag/blocks/spinhalf/compile.hpp>
#include <xdiag/blocks/spinhalf/spinhalf.hpp>
#include <xdiag/blocks/spinhalf/spinhalf_apply.hpp>
#include <xdiag/blocks/spinhalf/spinhalf_matrix.hpp>

#include <xdiag/blocks/electron/compile.hpp>
#include <xdiag/blocks/electron/electron.hpp>
#include <xdiag/blocks/electron/electron_apply.hpp>
#include <xdiag/blocks/electron/electron_matrix.hpp>

#include <xdiag/blocks/tj/compile.hpp>
#include <xdiag/blocks/tj/tj.hpp>
#include <xdiag/blocks/tj/tj_apply.hpp>
#include <xdiag/blocks/tj/tj_matrix.hpp>

#include <xdiag/states/fill.hpp>
#include <xdiag/states/gpwf.hpp>
#include <xdiag/states/product_state.hpp>
#include <xdiag/states/random_state.hpp>
#include <xdiag/states/state.hpp>

#include <xdiag/random/hash.hpp>
#include <xdiag/random/hash_functions.hpp>
#include <xdiag/random/random_utils.hpp>

#include <xdiag/combinatorics/binomial.hpp>
#include <xdiag/combinatorics/bit_patterns.hpp>
#include <xdiag/combinatorics/combinations.hpp>
#include <xdiag/combinatorics/combinations_index.hpp>
#include <xdiag/combinatorics/combinations_indexing.hpp>
#include <xdiag/combinatorics/fermi_table.hpp>
#include <xdiag/combinatorics/lin_table.hpp>
#include <xdiag/combinatorics/subsets.hpp>
#include <xdiag/combinatorics/subsets_index.hpp>
#include <xdiag/combinatorics/subsets_indexing.hpp>

#include <xdiag/basis/basis.hpp>

#include <xdiag/basis/spinhalf/basis_no_sz.hpp>
#include <xdiag/basis/spinhalf/basis_sublattice.hpp>
#include <xdiag/basis/spinhalf/basis_symmetric_no_sz.hpp>
#include <xdiag/basis/spinhalf/basis_symmetric_sz.hpp>
#include <xdiag/basis/spinhalf/basis_sz.hpp>

#include <xdiag/basis/tj/basis_np.hpp>
#include <xdiag/basis/tj/basis_symmetric_np.hpp>

#include <xdiag/basis/electron/basis_no_np.hpp>
#include <xdiag/basis/electron/basis_np.hpp>
#include <xdiag/basis/electron/basis_symmetric_no_np.hpp>
#include <xdiag/basis/electron/basis_symmetric_np.hpp>

#include <xdiag/symmetries/continuous_group.hpp>
#include <xdiag/symmetries/generated_group.hpp>
#include <xdiag/symmetries/group_action/group_action.hpp>
#include <xdiag/symmetries/group_action/group_action_lookup.hpp>
#include <xdiag/symmetries/group_action/group_action_sublattice.hpp>
#include <xdiag/symmetries/group_action/sublattice_stability.hpp>
#include <xdiag/symmetries/operations/fermi_sign.hpp>
#include <xdiag/symmetries/operations/group_action_operations.hpp>
#include <xdiag/symmetries/operations/representative_list.hpp>
#include <xdiag/symmetries/operations/symmetry_operations.hpp>
#include <xdiag/symmetries/permutation.hpp>
#include <xdiag/symmetries/permutation_group.hpp>
#include <xdiag/symmetries/qn.hpp>
#include <xdiag/symmetries/representation.hpp>

#include <xdiag/operators/bond.hpp>
#include <xdiag/operators/bondlist.hpp>
#include <xdiag/operators/bondlist_handler.hpp>
#include <xdiag/operators/compiler.hpp>
#include <xdiag/operators/non_branching_bonds.hpp>
#include <xdiag/operators/symmetrized_operator.hpp>


#include <xdiag/algorithms/norm_estimate.hpp>
#include <xdiag/algorithms/sparse_diag.hpp>

#include <xdiag/algorithms/lanczos/lanczos.hpp>
#include <xdiag/algorithms/lanczos/eigs_lanczos.hpp>
#include <xdiag/algorithms/lanczos/eigvals_lanczos.hpp>
#include <xdiag/algorithms/lanczos/lanczos_convergence.hpp>
#include <xdiag/algorithms/lanczos/tmatrix.hpp>
#include <xdiag/algorithms/lanczos/lanczos_pro.hpp>

#include <xdiag/algorithms/arnoldi/arnoldi.hpp>
#include <xdiag/algorithms/arnoldi/arnoldi_to_disk.hpp>
#include <xdiag/algorithms/gram_schmidt/gram_schmidt.hpp>
#include <xdiag/algorithms/gram_schmidt/orthogonalize.hpp>

#include <xdiag/algorithms/time_evolution/exp_sym_v.hpp>
#include <xdiag/algorithms/time_evolution/pade_matrix_exponential.hpp>
#include <xdiag/algorithms/time_evolution/time_evolution.hpp>
#include <xdiag/algorithms/time_evolution/zahexpv.hpp>

#include <xdiag/io/args.hpp>
#include <xdiag/io/file_h5.hpp>
#include <xdiag/io/file_toml.hpp>

#ifdef XDIAG_ENABLE_MPI
#include <xdiag/mpi/allreduce.hpp>
#include <xdiag/mpi/alltoall.hpp>
#include <xdiag/mpi/buffer.hpp>
#include <xdiag/mpi/datatype.hpp>
#include <xdiag/mpi/dot_mpi.hpp>
#include <xdiag/mpi/logger_mpi.hpp>
#include <xdiag/mpi/timing_mpi.hpp>
#include <xdiag/utils/print_mpi.hpp>

#include <xdiag/basis/spinhalf_mpi/spinhalf_mpi_basis_sz.hpp>
#include <xdiag/blocks/utils/block_utils_mpi.hpp>

#include <xdiag/blocks/spinhalf_mpi/spinhalf_mpi.hpp>
#include <xdiag/blocks/spinhalf_mpi/spinhalf_mpi_apply.hpp>
#include <xdiag/blocks/spinhalf_mpi/spinhalf_mpi_fill.hpp>
#include <xdiag/blocks/spinhalf_mpi/terms/get_prefix_postfix_mixed_bonds.hpp>

#include <xdiag/blocks/electron_mpi/electron_mpi.hpp>
#include <xdiag/blocks/electron_mpi/electron_mpi_apply.hpp>
#endif
