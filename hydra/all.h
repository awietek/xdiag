#pragma once

#include <lila/all.h>

#include "common.h"
#include "utils/bitops.h"
#include "utils/iochecks.h"

#include "combinatorics/binomial.h"
#include "combinatorics/subsets.h"
#include "combinatorics/bit_patterns.h"
#include "combinatorics/combinations.h"
#include "combinatorics/up_down_hole.h"

#include "indexing/lintable.h"

#include "models/models.h"
#include "models/utils/model_utils.h"
#include "models/utils/symmetrized_norm.h"

#include "models/spinhalf/spinhalf.h"
#include "models/spinhalf/spinhalf_matrix.h"
#include "models/spinhalf/spinhalf_apply.h"
#include "models/spinhalf/spinhalf_fill.h"

#include "models/spinhalf_symmetric/spinhalf_symmetric.h"
#include "models/spinhalf_symmetric/spinhalf_symmetric_matrix.h"
#include "models/spinhalf_symmetric/spinhalf_symmetric_apply.h"

#include "models/tj/tj.h"
#include "models/tj/tj_matrix.h"
#include "models/tj/tj_apply.h"

#include "models/tj_symmetric/tj_symmetric.h"
#include "models/tj_symmetric/tj_symmetric_matrix.h"
#include "models/tj_symmetric/tj_symmetric_apply.h"

#include "models/electron/electron.h"
#include "models/electron/electron_matrix.h"
#include "models/electron/electron_apply.h"

#include "models/electron_symmetric/electron_symmetric.h"
#include "models/electron_symmetric/electron_symmetric_matrix.h"
#include "models/electron_symmetric/electron_symmetric_apply.h"

#include "symmetries/permutation_group.h"
#include "symmetries/permutation_group_action.h"
#include "symmetries/permutation_group_lookup.h"
#include "symmetries/representation.h"
#include "symmetries/fermi_sign.h"
#include "symmetries/symmetry_utils.h"
#include "symmetries/symmetric_operator.h"

#include "operators/bond.h"
#include "operators/bondlist.h"
#include "operators/couplings.h"

#include "wavefunctions/gpwf_spinhalf.h"

#include "linalg/sparse_diag.h"
#include "linalg/algebra.h"
#include "linalg/lanczos/tmatrix.h"
#include "linalg/lanczos/lanczos_generic.h"
#include "linalg/lanczos/lanczos_convergence.h"
#include "linalg/lanczos/lanczos_eigenvalues.h"
#include "linalg/lanczos/lanczos_eigenvector.h"
