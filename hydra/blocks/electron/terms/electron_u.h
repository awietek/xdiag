#pragma once

#include <hydra/bitops/bitops.h>
#include <hydra/blocks/utils/block_utils.h>
#include <hydra/combinatorics/combinations.h>
#include <hydra/common.h>
#include <hydra/indexing/electron/electron_indexing.h>
#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>
#include <hydra/operators/operator_utils.h>

namespace hydra::terms {

template <typename bit_t, typename coeff_t, class Indexing, class Filler>
void electron_U(Couplings const &couplings, Indexing &&indexing,
                Filler &&fill) {
  using bitops::popcnt;

  if (couplings.defined("U")) {

    coeff_t U = utils::get_coupling<coeff_t>(couplings, "U");

    if (std::abs(U) > 1e-12) {
      idx_t idx = 0;
      for (bit_t up : indexing.states_ups()) {
        for (bit_t dn : indexing.states_dns()) {
          coeff_t val = U * (double)popcnt(up & dn);
          fill(idx, idx, val);
          ++idx;
        }
      }
    }
  }
}

} // namespace hydra::terms
