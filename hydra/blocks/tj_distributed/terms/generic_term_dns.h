#pragma once

#include <functional>
#include <hydra/bits/bitops.h>

namespace hydra::tj_distributed {

template <typename bit_t, typename coeff_t, bool fermi_ups,
          class BasisIn, class BasisOut, class NonZeroTermUps,
          class NonZeroTermDns, class TermAction, class Fill>
void generic_term_dns(BasisIn &&basis_in, BasisOut &&basis_out,
                      NonZeroTermUps &&non_zero_term_ups,
                      NonZeroTermDns &&non_zero_term_dns,
                      TermAction &&term_action, Fill &&fill) {

  (void)basis_in;
  (void)basis_out;
  (void)non_zero_term_ups;
  (void)non_zero_term_dns;
  (void)term_action;
  (void)fill;
}

} // namespace hydra::tj_distributed
