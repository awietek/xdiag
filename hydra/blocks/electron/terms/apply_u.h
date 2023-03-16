#pragma once

#include <hydra/bitops/bitops.h>
#include <hydra/common.h>

namespace hydra::electron {

template <typename bit_t, typename coeff_t, bool symmetric, class Indexing,
          class Filler>
void apply_u(coeff_t U, Indexing &&indexing, Filler &&fill) {

  if constexpr (symmetric) {
#ifdef _OPENMP
#pragma omp parallel for schedule(guided)
#endif
    for (idx_t idx_ups = 0; idx_ups < indexing.n_rep_ups(); ++idx_ups) {
      bit_t ups = indexing.rep_ups(idx_ups);
      auto dnss = indexing.dns_for_ups_rep(ups);
      idx_t idx = indexing.ups_offset(idx_ups);
      for (bit_t dns : dnss) {
        coeff_t val = U * (double)bitops::popcnt(ups & dns);
        fill(idx, idx, val);
        ++idx;
      }
    }
  } else {

#ifdef _OPENMP
#pragma omp parallel
    {
      auto ups_and_idces = indexing.states_indices_ups_thread();
#else
    auto ups_and_idces = indexing.states_indices_ups();
#endif

      for (auto [up, idx_up] : ups_and_idces) {

        idx_t idx = idx_up * indexing.size_dns();

        for (bit_t dn : indexing.states_dns()) {
          coeff_t val = U * (double)bitops::popcnt(up & dn);
          fill(idx, idx, val);
          ++idx;
        }
      }
#ifdef _OPENMP
    }
#endif
  }
}

} // namespace hydra::electron
