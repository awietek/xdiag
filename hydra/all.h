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

#ifndef HYDRA_HYDRA_H_
#define HYDRA_HYDRA_H_

#include "hilbertspaces/hubbard.h"
#include "hilbertspaces/spinhalf.h"

#include "indexing/indexhubbard.h"
#include "indexing/indextable.h"
#include "indexing/indexsearch.h"
#include "indexing/indexspinhalf.h"
#include "indexing/lintable.h"

#include "models/hubbardmodel.h"
#include "models/tjmodel.h"
#include "models/heisenbergmodel.h"
#include "models/spinlessfermions.h"

#include "utils/bitops.h"
#include "utils/combinatorics.h"
#include "utils/range.h"
#include "utils/typedefs.h"
#include "utils/iochecks.h"
#include "utils/format.h"
#include "utils/print.h"

#include "operators/bond.h"
#include "operators/bondlist.h"
#include "operators/couplings.h"

#include "thermodynamics/thermodynamics_detail.h"
#include "thermodynamics/thermodynamics_exact.h"
#include "thermodynamics/thermodynamics_tpq.h"

#include "dynamics/continuedfraction.h"

#include "symmetries/spacegroup.h"

#include "parameters/parser.h"
#include "parameters/parameters_impl.h"
#include "parameters/parameters.h"

namespace hydra 
{
  namespace all
  {
    using namespace hydra::dynamics;
    using namespace hydra::hilbertspaces;
    using namespace hydra::indexing;
    using namespace hydra::models;
    using namespace hydra::operators;
    using namespace hydra::parameters;
    using namespace hydra::symmetries;
    using namespace hydra::thermodynamics;
    using namespace hydra::utils;
    using namespace hydra::combinatorics;
  }
}


#endif
