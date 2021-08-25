#pragma once

#include <utility>

#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>
#include <hydra/symmetries/spacegroup.h>

namespace hydra {

template <class bit_t = std_bit_t,
          class SpaceGroupOperator = SpaceGroupOperator<bit_t>>
std::pair<BondList, Couplings>
SymmetricOperator(BondList const &bonds, Couplings const &cpls,
                  SpaceGroup const &space_group);

}
