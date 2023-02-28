#pragma once

#include <functional>

namespace hydra::electron {

template <typename bit_t, typename coeff_t, class Indexing, class TermAction,
          class Fill>
void apply_generic_term_diag_no_sym(Indexing &&indexing, TermAction &&term_action,
                            Fill &&fill) {
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
