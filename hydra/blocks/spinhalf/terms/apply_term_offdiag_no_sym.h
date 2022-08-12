#pragma once

#include <hydra/common.h>

namespace hydra::terms::spinhalf {

template <typename bit_t, typename coeff_t, class IndexingIn, class IndexingOut,
          class NonZeroTerm, class TermAction, class Fill>
void apply_term_offdiag_no_sym(IndexingIn &&indexing_in,
                               IndexingOut &&indexing_out,
                               NonZeroTerm &&non_zero_term,
                               TermAction &&term_action, Fill &&fill) {

  for (auto [spins_in, idx_in] : indexing_in) {

    if (non_zero_term(spins_in)) {
      auto [spins_out, coeff] = term_action(spins_in);
      auto idx_out = indexing_out.index(spins_out);
      assert(idx_out != invalid_index);
      fill(idx_out, idx_in, coeff);
    }
  }
}

} // namespace hydra::terms::spinhalf
