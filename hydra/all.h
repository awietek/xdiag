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
#include "extern/fmt/core.h"
#include "extern/clara/clara.hpp"

#include "common.h"
#include "config.h"

#include "bitops/bitops.h"
#include "utils/close.h"
#include "utils/iochecks.h"
#include "utils/logger.h"
#include "utils/precision.h"
#include "utils/print.h"
#include "utils/print_macro.h"
#include "utils/say_hello.h"

#include "blocks/blocks.h"
#include "blocks/utils/block_utils.h"

#include "algebra/algebra.h"
#include "algebra/matrix.h"

#include "states/gpwf_spinhalf.h"
#include "states/product_state.h"
#include "states/random_state.h"
#include "states/state.h"

#include "combinatorics/binomial.h"
#include "combinatorics/bit_patterns.h"
#include "combinatorics/combinations.h"
#include "combinatorics/combinations_index.h"
#include "combinatorics/subsets.h"
#include "combinatorics/subsets_index.h"

#include "indexing/combinations_indexing.h"
#include "indexing/fermi_table.h"
#include "indexing/indexing_variants.h"
#include "indexing/lin_table.h"
#include "indexing/subsets_indexing.h"

#include "indexing/spinhalf/indexing_no_sz.h"
#include "indexing/spinhalf/indexing_sublattice.h"
#include "indexing/spinhalf/indexing_symmetric_no_sz.h"
#include "indexing/spinhalf/indexing_symmetric_sz.h"
#include "indexing/spinhalf/indexing_sz.h"
#include "indexing/spinhalf/symmetric_iterator.h"

#include "indexing/tj/indexing_np.h"
#include "indexing/tj/indexing_symmetric_np.h"

#include "indexing/electron/indexing_no_np.h"
#include "indexing/electron/indexing_np.h"
#include "indexing/electron/indexing_symmetric_no_np.h"
#include "indexing/electron/indexing_symmetric_np.h"

#include "blocks/spinhalf/spinhalf.h"
#include "blocks/spinhalf/spinhalf_apply.h"
// #include "blocks/spinhalf/spinhalf_fill.h"
#include "blocks/spinhalf/spinhalf_matrix.h"
#include "blocks/spinhalf/terms/compile.h"

#include "blocks/electron/electron.h"
#include "blocks/electron/electron_apply.h"
#include "blocks/electron/electron_matrix.h"
#include "blocks/electron/terms/compile.h"

#include "blocks/tj/tj.h"
#include "blocks/tj/tj_apply.h"
#include "blocks/tj/tj_matrix.h"
#include "blocks/tj/tj_utils.h"

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
#include "symmetries/representation.h"

#include "operators/bond.h"
#include "operators/bondlist.h"
#include "operators/bondlist_handler.h"
#include "operators/compiler.h"
#include "operators/non_branching_bonds.h"
#include "operators/symmetrized_operator.h"

#include "algorithms/time_evolution/exp_sym_v.h"
#include "algorithms/lanczos/lanczos.h"
#include "algorithms/lanczos/lanczos_build.h"
#include "algorithms/lanczos/lanczos_convergence.h"
#include "algorithms/lanczos/lanczos_eigenvalues.h"
#include "algorithms/lanczos/lanczos_eigenvector.h"
#include "algorithms/lanczos/lanczos_vector_apply.h"
#include "algorithms/lanczos/tmatrix.h"
#include "algorithms/arnoldi/arnoldi.h"
#include "algorithms/arnoldi/arnoldi_to_disk.h"
#include "algorithms/sparse_diag.h"
#include "algorithms/time_evolution/time_evolution.h"
#include "algorithms/time_evolution/pade_matrix_exponential.h"
#include "algorithms/time_evolution/zahexpv.h"
#include "algorithms/norm_estimate.h"

#include "io/file_toml.h"
#include "io/file_h5.h"

#ifdef HYDRA_ENABLE_MPI
#include "mpi/allreduce.h"
#include "mpi/alltoall.h"
#include "mpi/buffer.h"
#include "mpi/datatype.h"
#include "mpi/dot_mpi.h"
#include "mpi/logger_mpi.h"
#include "mpi/timing_mpi.h"
#include "utils/print_mpi.h"

#include "blocks/utils/block_utils_mpi.h"

#include "indexing/spinhalf_mpi/spinhalf_mpi_indexing_sz.h"

#include "blocks/spinhalf_mpi/spinhalf_mpi.h"
#include "blocks/spinhalf_mpi/spinhalf_mpi_apply.h"
#include "blocks/spinhalf_mpi/spinhalf_mpi_fill.h"
#include "blocks/spinhalf_mpi/terms/get_prefix_postfix_mixed_bonds.h"

#include "blocks/electron_mpi/electron_mpi.h"
#include "blocks/electron_mpi/electron_mpi_apply.h"
#endif
