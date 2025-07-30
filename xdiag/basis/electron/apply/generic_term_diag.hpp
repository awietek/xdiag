// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/parallel/omp/omp_utils.hpp>

namespace xdiag::basis::electron {

template <typename bit_t, typename coeff_t, bool symmetric, class basis_t,
          class apply_f, class fill_f>
void generic_term_diag(Coupling cpl, basis_t const &basis, apply_f apply,
                       fill_f fill) try {
  if constexpr (symmetric) {
#ifdef _OPENMP
#pragma omp parallel
    {
      int num_thread = omp_get_thread_num();
#pragma omp for schedule(runtime)
#endif
      for (int64_t idx_ups = 0; idx_ups < basis.n_rep_ups(); ++idx_ups) {
        bit_t ups = basis.rep_ups(idx_ups);
        auto dnss = basis.dns_for_ups_rep(ups);
        int64_t idx = basis.ups_offset(idx_ups);
        for (bit_t dns : dnss) {
          coeff_t c = apply(ups, dns);
          XDIAG_FILL(idx, idx, c);
          ++idx;
        }
      }
#ifdef _OPENMP
    }
#endif

  } else { // if not symmetric
#ifdef _OPENMP
#pragma omp parallel
    {
      int num_thread = omp_get_thread_num();
      auto ups_and_idces = basis.states_indices_ups_thread();
#else
    auto ups_and_idces = basis.states_indices_ups();
#endif

      for (auto [ups, idx_up] : ups_and_idces) {
        int64_t idx = idx_up * basis.size_dns();
        for (bit_t dns : basis.states_dns()) {
          coeff_t c = apply(ups, dns);
          XDIAG_FILL(idx, idx, c);
          ++idx;
        }
      }
#ifdef _OPENMP
    }
#endif
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

} // namespace xdiag::basis::electron
