#pragma once

#include <lila/all.h>

#include "common.h"
#include "utils/bitops.h"
#include "utils/iochecks.h"

#include "combinatorics/binomial.h"
#include "combinatorics/bit_patterns.h"
#include "combinatorics/combinations.h"
#include "combinatorics/subsets.h"
#include "combinatorics/up_down_hole.h"

#include "indexing/indexing_symmetric.h"
#include "indexing/indexing_symmetric_fermionic.h"
#include "indexing/lintable.h"

#include "blocks/blocks.h"
#include "blocks/utils/block_utils.h"
#include "blocks/utils/symmetrized_norm.h"

#include "blocks/spinhalf/spinhalf.h"
#include "blocks/spinhalf/spinhalf_apply.h"
#include "blocks/spinhalf/spinhalf_fill.h"
#include "blocks/spinhalf/spinhalf_matrix.h"

#include "blocks/spinhalf_symmetric/spinhalf_symmetric.h"
#include "blocks/spinhalf_symmetric/spinhalf_symmetric_apply.h"
#include "blocks/spinhalf_symmetric/spinhalf_symmetric_matrix.h"

#include "blocks/tj/tj.h"
#include "blocks/tj/tj_apply.h"
#include "blocks/tj/tj_matrix.h"

// #include "blocks/tj_symmetric/tj_symmetric.h"
// #include "blocks/tj_symmetric/tj_symmetric_apply.h"
// #include "blocks/tj_symmetric/tj_symmetric_matrix.h"

#include "blocks/tj_symmetric_simple/tj_symmetric_simple.h"
#include "blocks/tj_symmetric_simple/tj_symmetric_simple_apply.h"
#include "blocks/tj_symmetric_simple/tj_symmetric_simple_matrix.h"

#include "blocks/electron/electron.h"
#include "blocks/electron/electron_apply.h"
#include "blocks/electron/electron_matrix.h"

#include "blocks/electron_symmetric/electron_symmetric.h"
#include "blocks/electron_symmetric/electron_symmetric_matrix.h"
#include "blocks/electron_symmetric/electron_symmetric_apply.h"

#include "blocks/electron_symmetric_simple/electron_symmetric_simple.h"
#include "blocks/electron_symmetric_simple/electron_symmetric_simple_apply.h"
#include "blocks/electron_symmetric_simple/electron_symmetric_simple_matrix.h"

#include "symmetries/fermi_sign.h"
#include "symmetries/permutation_group.h"
#include "symmetries/permutation_group_action.h"
#include "symmetries/permutation_group_lookup.h"
#include "symmetries/representation.h"
#include "symmetries/symmetric_operator.h"
#include "symmetries/symmetry_utils.h"

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
