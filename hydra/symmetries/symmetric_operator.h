#pragma once

#include <utility>

#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>
#include <hydra/symmetries/permutation_group.h>

namespace hydra {

std::pair<BondList, Couplings>
SymmetricOperator(BondList const &bonds, Couplings const &cpls,
                  PermutationGroup const & group);

}
