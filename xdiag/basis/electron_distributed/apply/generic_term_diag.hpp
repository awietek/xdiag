#pragma once

#include <xdiag/common.hpp>

namespace xdiag::basis::electron_distributed {

template <typename bit_t, typename coeff_t, class Basis, class TermAction>
void generic_term_diag(Basis &&basis, TermAction &&term_action,
                       const coeff_t *vec_in, coeff_t *vec_out) {
  int64_t idx = 0;
  for (bit_t up : basis.my_ups()) {
    for (bit_t dn : basis.all_dns()) {
      coeff_t val = term_action(up, dn);
      vec_out[idx] += val * vec_in[idx];
      ++idx;
    }
  }
}

} // namespace xdiag::basis::electron_distributed
