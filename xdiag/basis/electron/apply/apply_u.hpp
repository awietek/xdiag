#pragma once

#include <xdiag/bits/bitops.hpp>
#include <xdiag/common.hpp>

namespace xdiag::basis::electron {

template <typename bit_t, typename coeff_t, bool symmetric, class basis_t,
          class fill_f>
void apply_u(Coupling const &cpl, basis_t const &basis, fill_f fill) try {

  coeff_t U = cpl.scalar().as<coeff_t>();

  if constexpr (symmetric) {
#ifdef _OPENMP
#pragma omp parallel for schedule(guided)
#endif
    for (int64_t idx_ups = 0; idx_ups < basis.n_rep_ups(); ++idx_ups) {
      bit_t ups = basis.rep_ups(idx_ups);
      auto dnss = basis.dns_for_ups_rep(ups);
      int64_t idx = basis.ups_offset(idx_ups);
      for (bit_t dns : dnss) {
        coeff_t val = U * (double)bits::popcnt(ups & dns);
        fill(idx, idx, val);
        ++idx;
      }
    }
  } else {

#ifdef _OPENMP
#pragma omp parallel
    {
      auto ups_and_idces = basis.states_indices_ups_thread();
#else
    auto ups_and_idces = basis.states_indices_ups();
#endif

      for (auto [up, idx_up] : ups_and_idces) {
        int64_t idx = idx_up * basis.size_dns();
        for (bit_t dn : basis.states_dns()) {
          coeff_t val = U * (double)bits::popcnt(up & dn);
          fill(idx, idx, val);
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
