#pragma once

#include <vector>

namespace hydra::tj {

template <typename bit_t, typename coeff_t, bool symmetric, class Indexing,
          class TermAction, class Fill>
void generic_term_diag(IndexingIn &&indexing, TermAction &&term_action,
                       Fill &&fill) {
  int n_sites = indexing.n_sites();
  bit_t sitesmask = ((bit_t)1 << n_sites) - 1;

  if constexpr (symmetric) {
    idx_t idx = 0;
    for (idx_t idx_up_in = 0; idx_up_in < indexing.n_rep_ups(); ++idx_up_in) {
      bit_t up_in = indexing.rep_ups(idx_up_in);
      bit_t not_up_in = (~up_in) & sitesmask;

      auto dncs_in = indexing_in.dns_for_ups_rep(up_in);
      for (bit_t dnc_in : dncs_in) {
        idx_t dn_in = bitops::deposit(dnc_in, not_up_in);
        coeff_t val = term_action(up_in, dn_in);
        fill(idx, idx, val);
        ++idx
      }
    }
  }

  else { // if not symmetric
    auto ups_and_idces = indexing_in.states_indices_ups();
    for (auto [up_in, idx_up_in] : ups_and_idces) {
      bit_t not_up_in = (~up_in) & sitesmask;
      auto dncs_in = indexing_in.states_dncs(up_in);
      for (bit_t dnc_in : dncs_in) {
        bit_t dn_in = bitops::deposit(dnc_in, not_up_in);
        coeff_t val = term_action(up_in, dn_in);
	++idx;
      }
    }
  }
} // namespace hydra::tj
