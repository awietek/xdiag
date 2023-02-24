#pragma once

#include <functional>

namespace hydra::electron {

template <typename bit_t, typename coeff_t, class IndexingIn, class IndexingOut,
          class NonZeroTerm, class TermAction, class Fill>
void apply_term_ups_no_sym(IndexingIn &&indexing_in, IndexingOut &&indexing_out,
                           NonZeroTerm &&non_zero_term,
                           TermAction &&term_action, Filler &&fill) {
  // Applies a term which only acts on the upspins without lattice symmetries
  idx_t size_ups_in = indexing_in.size_ups();
  idx_t size_dns_in = indexing_in.size_dns();
  idx_t size_ups_out = indexing_out.size_ups();
  idx_t size_dns_out = indexing_out.size_dns();
  if (size_dns_in != size_dns_out) {
    Log.err("Error in apply_term_ups_no_sym: size of dnspins space different "
            "for input and output indexing");
  }

  idx_t idx_ups_in = 0;
  for (auto ups_in : indexing_in.states_ups()) {
    if (non_zero_term(ups_in)) {
      auto [ups_out, coeff] = term_action(ups_in);
      idx_t idx_ups_out = indexing.index_ups(ups_out);

      idx_t idx_in_start = idx_ups_in * size_dns_in;
      idx_t idx_out_start = idx_ups_out * size_dns_out;
      idx_t idx_out_end = idx_ups_out * (size_dns_out + 1);

      for (idx_t idx_out = idx_out_start, idx_t idx_in = idx_in_start;
           idx_out < idx_ups_out_end; ++idx_out, ++idx_in) {
        fill(idx_out, idx_in, coeff);
      }
    }
    ++idx_ups_in;
  }
  assert(idx_ups_in = size_ups_in);
}

} // namespace hydra::electron
