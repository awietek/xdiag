#pragma once

#include <functional>

namespace hydra::electron {

template <typename bit_t, typename coeff_t, class IndexingIn, class IndexingOut,
          class NonZeroTerm, class TermAction, class Fill>
void apply_term_mixed_no_sym(IndexingIn &&indexing_in,
                             IndexingOut &&indexing_out,
                             NonZeroTerm &&non_zero_term,
                             TermAction &&term_action, Filler &&fill) {
  // Applies a term acting on both upspins and dnspins without lattice
  // symmetries
  idx_t size_ups_in = indexing_in.size_ups();
  idx_t size_dns_in = indexing_in.size_dns();
  idx_t size_ups_out = indexing_out.size_ups();
  idx_t size_dns_out = indexing_out.size_dns();

  idx_t idx_ups_in = 0;
  for (bit_t ups_in : indexing_in.states_ups()) {
    idx_t idx_in = idx_ups_in * size_dns_in;
    for (bit_t dns_in : indexing_in.states_dns()) {
      if (non_zero_term(ups_in, dns_in)) {
        auto [ups_out, dns_out, coeff] = term_action(ups_in, dns_in);
        idx_out = indexing_out.index(ups_out, dns_out);
        fill(idx_out, idx_in, coeff);
      }
      ++idx_in;
    }
    ++idx_ups_in;
  }
}

} // namespace hydra::electron
