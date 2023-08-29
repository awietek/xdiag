#pragma once

#include <vector>
#include <hydra/bits/bitops.h>

namespace hydra::tj {

template <typename bit_t, typename coeff_t, bool symmetric, class Basis,
          class TermAction, class Fill>
void generic_term_diag(Basis &&basis, TermAction &&term_action,
                       Fill &&fill) {
  int64_t n_sites = basis.n_sites();
  bit_t sitesmask = ((bit_t)1 << n_sites) - 1;

  if constexpr (symmetric) {

#ifdef _OPENMP
#pragma omp parallel for schedule(guided)
#endif
    for (idx_t idx_up_in = 0; idx_up_in < basis.n_rep_ups(); ++idx_up_in) {
      idx_t idx = basis.ups_offset(idx_up_in);

      bit_t up_in = basis.rep_ups(idx_up_in);
      bit_t not_up_in = (~up_in) & sitesmask;

      auto dnss_in = basis.dns_for_ups_rep(up_in);
      auto up_syms = basis.syms_ups(up_in);

      // Trivial stabilizer of ups -> dns need to be deposited
      if (up_syms.size() == 1) {
        for (bit_t dnc_in : dnss_in) {
          idx_t dn_in = bits::deposit(dnc_in, not_up_in);
          coeff_t val = term_action(up_in, dn_in);
          fill(idx, idx, val);
          ++idx;
        }
      }

      // stabilizer for ups -> dns don't need to be deposited
      else {
        for (bit_t dn_in : dnss_in) {
          coeff_t val = term_action(up_in, dn_in);
          fill(idx, idx, val);
          ++idx;
        }
      }
    }
  }

  else { // if not symmetric
#ifdef _OPENMP
#pragma omp parallel
    {
      auto ups_and_idces = basis.states_indices_ups_thread();
#else
    auto ups_and_idces = basis.states_indices_ups();
#endif
      for (auto [up_in, idx_up_in] : ups_and_idces) {
        bit_t not_up_in = (~up_in) & sitesmask;
        auto dncs_in = basis.states_dncs(up_in);
        idx_t idx = basis.ups_offset(idx_up_in);
        for (bit_t dnc_in : dncs_in) {
          bit_t dn_in = bits::deposit(dnc_in, not_up_in);
          coeff_t val = term_action(up_in, dn_in);
          fill(idx, idx, val);
          ++idx;
        }
      }

#ifdef _OPENMP
    }
#endif
  } // if not symmetric
}
} // namespace hydra::tj
