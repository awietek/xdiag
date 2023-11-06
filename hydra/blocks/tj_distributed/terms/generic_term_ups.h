#pragma once

#include <hydra/bits/bitops.h>
#include <vector>

namespace hydra::tj_distributed {

template <typename bit_t, typename coeff_t, class BasisIn, class BasisOut,
          class NonZeroTerm, class TermAction, class Fill>
void generic_term_ups(BasisIn &&basis_in, BasisOut &&basis_out,
                      NonZeroTerm &&non_zero_term, TermAction &&term_action,
                      Fill &&fill) {
  (void)basis_in;
  (void)basis_out;
  (void)non_zero_term;
  (void)term_action;
  (void)fill;
}

} // namespace hydra::tj_distributed
