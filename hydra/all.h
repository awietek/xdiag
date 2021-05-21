// Copyright 2019 Alexander Wietek - All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#pragma once

#include "common.h"
#include "utils/bitops.h"
#include "utils/iochecks.h"

#include "combinatorics/binomial.h"
#include "combinatorics/subsets.h"
#include "combinatorics/bit_patterns.h"
#include "combinatorics/combinations.h"
#include "combinatorics/up_down_hole.h"

#include "indexing/lintable.h"

#include "models/spinhalf/spinhalf.h"

#include "models/electron/electron.h"
#include "models/electron/electron_matrix.h"
#include "models/electron/electron_apply.h"
#include "models/electron/electron_symmetric.h"
#include "models/electron/electron_symmetric_matrix.h"
#include "models/electron/electron_symmetric_apply.h"

#include "models/tj/tj.h"
#include "models/tj/tj_matrix.h"
#include "models/tj/tj_apply.h"

#include "symmetries/spacegroup.h"
#include "symmetries/spinflip.h"
#include "symmetries/representation.h"
#include "symmetries/spacegroup_operations.h"
#include "symmetries/spacegroup_operator.h"
#include "symmetries/fermi_sign.h"

#include "operators/bond.h"
#include "operators/bondlist.h"
#include "operators/couplings.h"

#include "algebra/diagonalization.h"

// #include "combinatorics/up_down_hole.h"
// #include "symmetries/charactertable.h"


// #include "qns/qn_spinhalf.h"
// #include "qns/qn_electron.h"
// #include "qns/qn_tj.h"

// #include "states/state_spinhalf.h"
// #include "states/state_electron.h"
// #include "states/state_tj.h"

// #include "bases/basis_spinhalf.h"
// #include "bases/basis_tj.h"
// #include "bases/basis_electron.h"

// #include "indexing/index_table.h"
// #include "indexing/index_search.h"
// #include "indexing/index_spinhalf.h"
// #include "indexing/index_electron.h"

// #include "models/hubbardmodel.h"
// #include "models/tjmodel.h"
// #include "models/heisenbergmodel.h"
// #include "models/spinlessfermions.h"



// #include "utils/range.h"
// #include "utils/iochecks.h"
// #include "utils/format.h"
// #include "utils/print.h"
// #include "utils/hash.h"



// #include "entanglement/reduced_density_matrix.h"
// #include "entanglement/entanglement_entropy.h"

// // #include "thermodynamics/thermodynamics_detail.h"
// // #include "thermodynamics/thermodynamics_exact.h"
// // #include "thermodynamics/thermodynamics_tpq.h"
// // #include "dynamics/continuedfraction.h"

// 

// #include "parameters/parser.h"
// #include "parameters/parameters_impl.h"
// #include "parameters/parameters.h"
