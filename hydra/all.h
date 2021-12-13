#pragma once

#include <lila/all.h>

#include "bitops/bitops.h"
#include "common.h"
#include "utils/iochecks.h"
#include "utils/openmp_utils.h"

#include "combinatorics/binomial.h"
#include "combinatorics/bit_patterns.h"
#include "combinatorics/combinations.h"
#include "combinatorics/combinations_index.h"
#include "combinatorics/subsets.h"
#include "combinatorics/subsets_index.h"

#include "indexing/lintable.h"
#include "indexing/lintable_raw.h"
#include "indexing/spinhalf/spinhalf_indexing.h"
#include "indexing/spinhalf/spinhalf_indexing_no_sz.h"
#include "indexing/spinhalf/spinhalf_symmetric_indexing.h"
#include "indexing/tj/tj_indexing.h"
#include "indexing/tj/tj_symmetric_indexing.h"

#include "blocks/blocks.h"
#include "blocks/utils/block_utils.h"

#include "blocks/spinhalf/spinhalf.h"
#include "blocks/spinhalf/spinhalf_apply.h"
#include "blocks/spinhalf/spinhalf_fill.h"
#include "blocks/spinhalf/spinhalf_matrix.h"

#include "blocks/tj/tj.h"
#include "blocks/tj/tj_apply.h"
#include "blocks/tj/tj_matrix.h"
#include "blocks/tj/tj_utils.h"

#include "blocks/tj_symmetric/tj_symmetric.h"
#include "blocks/tj_symmetric/tj_symmetric_apply.h"
#include "blocks/tj_symmetric/tj_symmetric_matrix.h"

#include "blocks/electron/electron.h"
#include "blocks/electron/electron_apply.h"
#include "blocks/electron/electron_matrix.h"

#include "blocks/electron_symmetric/electron_symmetric.h"
#include "blocks/electron_symmetric/electron_symmetric_apply.h"
#include "blocks/electron_symmetric/electron_symmetric_matrix.h"

#include "symmetries/fermi_sign.h"
#include "symmetries/permutation_group.h"
#include "symmetries/permutation_group_action.h"
#include "symmetries/permutation_group_lookup.h"
#include "symmetries/representation.h"
#include "symmetries/symmetric_operator.h"
#include "symmetries/symmetry_operations.h"
#include "symmetries/representative_list.h"
#include "symmetries/fermi_bool_table.h"

#include "operators/bond.h"
#include "operators/bondlist.h"
#include "operators/couplings.h"

#include "wavefunctions/gpwf_spinhalf.h"

#include "linalg/algebra.h"
#include "linalg/lanczos/lanczos_convergence.h"
#include "linalg/lanczos/lanczos_eigenvalues.h"
#include "linalg/lanczos/lanczos_eigenvector.h"
#include "linalg/lanczos/lanczos_generic.h"
#include "linalg/lanczos/tmatrix.h"
#include "linalg/sparse_diag.h"
