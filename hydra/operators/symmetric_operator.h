#pragma once

#include <utility>

#include <hydra/operators/bondlist.h>
#include <hydra/symmetries/permutation_group.h>

namespace hydra {

BondList symmetric_operator(BondList const &bonds,
                            PermutationGroup const &group);

}
