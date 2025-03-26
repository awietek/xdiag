#pragma once

#include <xdiag/common.hpp>

namespace xdiag::basis::tj_distributed {

template <typename bit_t, typename coeff_t, class basis_t, class term_action_f>
void generic_term_diag(basis_t const &basis, term_action_f term_action,
                       const coeff_t *vec_in, coeff_t *vec_out) {
  int64_t idx_up = 0;
  int64_t idx = 0;
  for (bit_t up : basis.my_ups()) {
    for (bit_t dn : basis.my_dns_for_ups(idx_up)) {
      coeff_t val = term_action(up, dn);
      vec_out[idx] += val * vec_in[idx];
      ++idx;
    }
    ++idx_up;
  }
}

} // namespace xdiag::basis::tj_distributed
