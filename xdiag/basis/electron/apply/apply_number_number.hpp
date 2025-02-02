#pragma once

#include <xdiag/bits/gbit.hpp>

namespace xdiag::basis::electron {

template <typename bit_t, typename coeff_t, bool symmetric, class Basis,
          class Fill>
void apply_number_number(Coupling const &cpl, Op const &op, Basis &&basis,
                         Fill fill) try {
  using bits::gbit;

  coeff_t mu = cpl.scalar().as<coeff_t>();
  int64_t s1 = op[0];
  int64_t s2 = op[1];

  if constexpr (symmetric) {
    int64_t idx = 0;
    for (int64_t idx_ups = 0; idx_ups < basis.n_rep_ups(); ++idx_ups) {
      bit_t ups = basis.rep_ups(idx_ups);
      auto dnss = basis.dns_for_ups_rep(ups);
      for (bit_t dns : dnss) {
        int n1 = gbit(ups, s1) + gbit(dns, s1);
        int n2 = gbit(ups, s2) + gbit(dns, s2);
        fill(idx, idx, mu * (coeff_t)(n1 * n2));
        ++idx;
      }
    }

  } else { // if not symmetric
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
          int n1 = gbit(up, s1) + gbit(dn, s1);
          int n2 = gbit(up, s2) + gbit(dn, s2);
          fill(idx, idx, mu * (coeff_t)(n1 * n2));
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
