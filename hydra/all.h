#pragma once

#if __cplusplus < 201703L
#error "Hydra requires at least C++17"
#endif

#ifdef HYDRA_USE_MPI
#include <mpi.h>
#endif

#ifdef _OPENMP
#include "parallel/omp/omp_utils.h"
#include <omp.h>
// #include "symmetries/operations/representative_list_omp.h"
#endif

#ifdef HYDRA_USE_HDF5
#define ARMA_USE_HDF5
#endif

#include "extern/armadillo/armadillo"
#include "extern/clara/clara.hpp"
#include "extern/fmt/core.h"

#include "common.h"
#include "config.h"

#include "bits/bitops.h"
#include "utils/close.h"
#include "utils/iochecks.h"
#include "utils/logger.h"
#include "utils/precision.h"
#include "utils/print.h"
#include "utils/print_macro.h"
#include "utils/say_hello.h"

#include "algebra/algebra.h"
#include "algebra/apply.h"
#include "algebra/matrix.h"

// Includes for different block types
#include "blocks/blocks.h"

#include "blocks/spinhalf/compile.h"
#include "blocks/spinhalf/spinhalf.h"
#include "blocks/spinhalf/spinhalf_apply.h"
#include "blocks/spinhalf/spinhalf_matrix.h"

#include "blocks/electron/compile.h"
#include "blocks/electron/electron.h"
#include "blocks/electron/electron_apply.h"
#include "blocks/electron/electron_matrix.h"

#include "blocks/tj/compile.h"
#include "blocks/tj/tj.h"
#include "blocks/tj/tj_apply.h"
#include "blocks/tj/tj_matrix.h"

#include "states/fill.h"
#include "states/gpwf.h"
#include "states/product_state.h"
#include "states/random_state.h"
#include "states/state.h"

#include "random/hash.h"
#include "random/hash_functions.h"
#include "random/random_utils.h"

#include "combinatorics/binomial.h"
#include "combinatorics/bit_patterns.h"
#include "combinatorics/combinations.h"
#include "combinatorics/combinations_index.h"
#include "combinatorics/combinations_indexing.h"
#include "combinatorics/fermi_table.h"
#include "combinatorics/lin_table.h"
#include "combinatorics/subsets.h"
#include "combinatorics/subsets_index.h"
#include "combinatorics/subsets_indexing.h"

#include "basis/basis.h"

#include "basis/spinhalf/basis_no_sz.h"
#include "basis/spinhalf/basis_sublattice.h"
#include "basis/spinhalf/basis_symmetric_no_sz.h"
#include "basis/spinhalf/basis_symmetric_sz.h"
#include "basis/spinhalf/basis_sz.h"

#include "basis/tj/basis_np.h"
#include "basis/tj/basis_symmetric_np.h"

#include "basis/electron/basis_no_np.h"
#include "basis/electron/basis_np.h"
#include "basis/electron/basis_symmetric_no_np.h"
#include "basis/electron/basis_symmetric_np.h"

#include "symmetries/continuous_group.h"
#include "symmetries/generated_group.h"
#include "symmetries/group_action/group_action.h"
#include "symmetries/group_action/group_action_lookup.h"
#include "symmetries/group_action/group_action_sublattice.h"
#include "symmetries/group_action/sublattice_stability.h"
#include "symmetries/operations/fermi_sign.h"
#include "symmetries/operations/group_action_operations.h"
#include "symmetries/operations/representative_list.h"
#include "symmetries/operations/symmetry_operations.h"
#include "symmetries/permutation.h"
#include "symmetries/permutation_group.h"
#include "symmetries/qn.h"
#include "symmetries/representation.h"

#include "operators/bond.h"
#include "operators/bondlist.h"
#include "operators/bondlist_handler.h"
#include "operators/compiler.h"
#include "operators/non_branching_bonds.h"
#include "operators/symmetrized_operator.h"

#include "algorithms/lanczos/eigs_lanczos.h"
#include "algorithms/lanczos/eigvals_lanczos.h"
#include "algorithms/lanczos/lanczos_convergence.h"
#include "algorithms/lanczos/tmatrix.h"
#include "algorithms/sparse_diag.h"

// #include "algorithms/arnoldi/arnoldi.h"
// #include "algorithms/arnoldi/arnoldi_to_disk.h"
// #include "algorithms/gram_schmidt/gram_schmidt.h"
// #include "algorithms/gram_schmidt/orthogonalize.h"
// #include "algorithms/lanczos/lanczos.h"
// #include "algorithms/lanczos/lanczos_build.h"
// #include "algorithms/lanczos/lanczos_convergence.h"
// #include "algorithms/lanczos/lanczos_eigenvalues.h"
// #include "algorithms/lanczos/lanczos_eigenvector.h"
// #include "algorithms/lanczos/lanczos_pro.h"
// #include "algorithms/lanczos/lanczos_vector_apply.h"
// #include "algorithms/lanczos/tmatrix.h"
// #include "algorithms/norm_estimate.h"
// #include "algorithms/sparse_diag.h"
// #include "algorithms/time_evolution/exp_sym_v.h"
// #include "algorithms/time_evolution/pade_matrix_exponential.h"
// #include "algorithms/time_evolution/time_evolution.h"
// #include "algorithms/time_evolution/zahexpv.h"

#include "io/args.h"
#include "io/file_h5.h"
#include "io/file_toml.h"

#ifdef HYDRA_ENABLE_MPI
#include "mpi/allreduce.h"
#include "mpi/alltoall.h"
#include "mpi/buffer.h"
#include "mpi/datatype.h"
#include "mpi/dot_mpi.h"
#include "mpi/logger_mpi.h"
#include "mpi/timing_mpi.h"
#include "utils/print_mpi.h"

#include "basis/spinhalf_mpi/spinhalf_mpi_basis_sz.h"
#include "blocks/utils/block_utils_mpi.h"

#include "blocks/spinhalf_mpi/spinhalf_mpi.h"
#include "blocks/spinhalf_mpi/spinhalf_mpi_apply.h"
#include "blocks/spinhalf_mpi/spinhalf_mpi_fill.h"
#include "blocks/spinhalf_mpi/terms/get_prefix_postfix_mixed_bonds.h"

#include "blocks/electron_mpi/electron_mpi.h"
#include "blocks/electron_mpi/electron_mpi_apply.h"
#endif
