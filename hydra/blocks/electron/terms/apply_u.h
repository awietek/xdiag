#pragma once

#include <hydra/bitops/bitops.h>
#include <hydra/common.h>

namespace hydra::electron {

template <typename bit_t, typename coeff_t, class Indexing, class Filler>
void apply_u(coeff_t U, Indexing &&indexing, Filler &&fill) {
  idx_t idx = 0;
  for (bit_t up : indexing.states_ups()) {
    for (bit_t dn : indexing.states_dns()) {
      coeff_t val = U * (double)bitops::popcnt(up & dn);
      fill(idx, idx, val);
      ++idx;
    }
  }
}

} // namespace hydra::electron
