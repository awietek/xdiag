#pragma once

#ifdef HYDRA_ENABLE_MPI
#include <mpi.h>
#endif

#ifdef HYDRA_ENABLE_OPENMP
#include <omp.h>
#endif

#include <lila/all.h>

#include "bitops/bitops.h"
#include "common.h"
#include "utils/iochecks.h"
#include "utils/openmp_utils.h"
#include "utils/print.h"
#include "utils/logger.h"
#include "utils/print_macro.h"

#include "combinatorics/binomial.h"
#include "combinatorics/bit_patterns.h"
#include "combinatorics/combinations.h"
#include "combinatorics/combinations_index.h"
#include "combinatorics/combinatorics_omp_utils.h"
#include "combinatorics/subsets.h"
#include "combinatorics/subsets_index.h"

#include "indexing/combinations_indexing.h"
#include "indexing/lintable.h"
#include "indexing/subsets_indexing.h"
#include "indexing/indexing_variants.h"


#include "indexing/spinhalf/spinhalf_indexing_sz.h"
#include "indexing/spinhalf/spinhalf_indexing_no_sz.h"
#include "indexing/spinhalf/spinhalf_symmetric_indexing_sz.h"
#include "indexing/spinhalf/spinhalf_symmetric_indexing_no_sz.h"

#include "indexing/tj/tj_indexing.h"
#include "indexing/tj/tj_symmetric_indexing.h"

#include "indexing/electron/electron_indexing.h"
#include "indexing/electron/electron_indexing_no_np.h"
#include "indexing/electron/electron_symmetric_indexing.h"
#include "indexing/electron/electron_symmetric_indexing_no_np.h"

#include "blocks/blocks.h"
#include "blocks/target_block.h"
#include "blocks/utils/block_utils.h"

#include "blocks/spinhalf/spinhalf.h"
#include "blocks/spinhalf/spinhalf_apply.h"
#include "blocks/spinhalf/spinhalf_fill.h"
#include "blocks/spinhalf/spinhalf_matrix.h"

#include "blocks/electron/electron.h"
#include "blocks/electron/electron_apply.h"
#include "blocks/electron/electron_matrix.h"

#include "blocks/tj/tj.h"
#include "blocks/tj/tj_apply.h"
#include "blocks/tj/tj_matrix.h"
#include "blocks/tj/tj_utils.h"

#include "symmetries/fermi_bool_table.h"
#include "symmetries/fermi_sign.h"
#include "symmetries/permutation_group.h"
#include "symmetries/representation.h"
#include "symmetries/representative_list.h"
#include "symmetries/symmetric_operator.h"
#include "symmetries/symmetry_operations.h"
#include "symmetries/group_action/group_action.h"
#include "symmetries/group_action/group_action_lookup.h"

#include "operators/bond.h"
#include "operators/bondlist.h"
#include "operators/couplings.h"
#include "operators/operator_utils.h"
#include "operators/operator_qns.h"

#include "wavefunctions/gpwf_spinhalf.h"

#include "linalg/algebra.h"
#include "linalg/state.h"
#include "linalg/lanczos/lanczos_convergence.h"
#include "linalg/lanczos/lanczos_eigenvalues.h"
#include "linalg/lanczos/lanczos_eigenvector.h"
#include "linalg/lanczos/lanczos_generic.h"
#include "linalg/lanczos/tmatrix.h"
#include "linalg/sparse_diag.h"

#ifdef HYDRA_ENABLE_MPI
#include "utils/print_mpi.h"
#include "mpi/datatype.h"
#include "mpi/allreduce.h"
#include "mpi/alltoall.h"
#include "mpi/logger_mpi.h"
#include "mpi/timing_mpi.h"
#include "mpi/dot_mpi.h"
#include "mpi/buffer.h"

#include "blocks/utils/block_utils_mpi.h"

#include "indexing/spinhalf_mpi/spinhalf_mpi_indexing_sz.h"

#include "blocks/spinhalf_mpi/terms/get_prefix_postfix_mixed_bonds.h"
#include "blocks/spinhalf_mpi/spinhalf_mpi.h"
#include "blocks/spinhalf_mpi/spinhalf_mpi_apply.h"
#include "blocks/spinhalf_mpi/spinhalf_mpi_fill.h"

#include "blocks/electron_mpi/electron_mpi.h"
#include "blocks/electron_mpi/electron_mpi_apply.h"
#endif
