#pragma once

#include <xdiag/common.hpp>

namespace xdiag::basis::electron_distributed {

template <typename coeff_t, class basis_t, class term_action_f>
void generic_term_diag(basis_t const &basis, term_action_f term_action,
                       const coeff_t *vec_in, coeff_t *vec_out) {
  using bit_t = typename basis_t::bit_t;

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
