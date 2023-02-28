#pragma once

#include <functional>

namespace hydra::electron {

template <typename bit_t, typename coeff_t, bool fermi_ups, class IndexingIn,
          class IndexingOut, class NonZeroTerm, class TermAction, class Fill>
void apply_generic_term_dns_no_sym(IndexingIn &&indexing_in,
                                   IndexingOut &&indexing_out,
                                   NonZeroTerm &&non_zero_term,
                                   TermAction &&term_action, Fill &&fill) {
  // Applies a term which only acts on the dnspins without lattice symmetries
  idx_t size_ups_in = indexing_in.size_ups();
  idx_t size_ups_out = indexing_out.size_ups();

  idx_t size_dns_in = indexing_in.size_dns();
  idx_t size_dns_out = indexing_out.size_dns();
  if (size_ups_in != size_ups_out) {
    Log.err("Error in apply_term_dns_no_sym: size of upspins space different "
            "for input and output indexing");
  }

  idx_t idx_dns_in = 0;
  for (bit_t dns_in : indexing_in.states_dns()) {
    if (non_zero_term(dns_in)) {
      auto [dns_out, coeff] = term_action(dns_in);

      idx_t idx_dns_out = indexing_out.index_dns(dns_out);
      idx_t idx_out_start = idx_dns_out;
      idx_t idx_out_end = idx_dns_out + indexing_out.size();

      if constexpr (fermi_ups) {
        idx_t idx_out = idx_out_start;
        idx_t idx_in = idx_dns_in;
        for (bit_t ups : indexing_in.states_ups()) {
          bool fermi = bitops::popcnt(ups) & 1;
          fill(idx_out, idx_in, fermi ? coeff : -coeff);
          idx_out += size_dns_out;
          idx_in += size_dns_in;
        }
      } else {
        for (idx_t idx_out = idx_out_start, idx_in = idx_dns_in;
             idx_out < idx_out_end;
             idx_out += size_dns_out, idx_in += size_dns_in) {
          fill(idx_out, idx_in, coeff);
        }
      }
    }
    ++idx_dns_in;
  }
  assert(idx_dns_in = size_dns_in);
}

} // namespace hydra::electron
