#pragma once

#include <functional>

namespace hydra::electron {

template <typename bit_t, typename coeff_t, class Indexing, class TermAction,
          class Fill>
void apply_term_diag(Indexing &&indexing, NonZeroTerm &&non_zero_term,
                     TermAction &&term_action, Filler &&fill) {
  idx_t idx = 0;
  for (auto ups : indexing.states_ups()) {
    for (auto dns : indexing.states_dns()) {
      coeff_t coeff = term_action(ups, dns);
      fill(idx, idx, coeff);
      ++idx;
    }
  }
}

} // namespace hydra::electron
