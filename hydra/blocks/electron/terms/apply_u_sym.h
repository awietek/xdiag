#pragma once

#include <hydra/bitops/bitops.h>
#include <hydra/common.h>

namespace hydra::electron {

template <typename bit_t, typename coeff_t, class Indexing, class Filler>
void apply_u_sym(coeff_t U, Indexing &&indexing, Filler &&fill) {

  idx_t idx = 0;
  // #ifdef _OPENMP
  // #pragma omp parallel for schedule(guided)
  // #endif
  for (idx_t idx_ups = 0; idx_ups < indexing.n_rep_ups(); ++idx_ups) {
    bit_t ups = indexing.rep_ups(idx_ups);
    auto dnss = indexing.dns_for_ups_rep(ups);
    for (bit_t dns : dnss) {
      coeff_t val = U * (double)bitops::popcnt(ups & dns);
      fill(idx, idx, val);
      ++idx;
    }
  }
}

} // namespace hydra::electron
