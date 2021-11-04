#pragma once

#include <lila/utils/logger.h>

#include <hydra/bitops/bitops.h>
#include <hydra/blocks/tj/tj.h>
#include <hydra/blocks/utils/block_utils.h>
#include <hydra/combinatorics/combinations.h>
#include <hydra/common.h>
#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>

namespace hydra::terms::tj {


template <class bit_t, class Filler>
void do_exchange(BondList const &bonds, Couplings const &couplings,
                 tJ<bit_t> const &block, Filler &&fill) {
}

} // namespace hydra::terms::tj
