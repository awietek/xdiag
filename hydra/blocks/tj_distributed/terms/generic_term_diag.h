#pragma once

#include <hydra/common.h>

namespace hydra::tj_distributed {

template <typename bit_t, typename coeff_t, class Basis, class TermAction>
void generic_term_diag(Basis &&basis, TermAction &&term_action,
                       arma::Col<coeff_t> const &vec_in,
                       arma::Col<coeff_t> &vec_out) {
  int64_t idx_up = 0;
  int64_t idx = 0;
  for (bit_t up : basis.my_ups()) {
    for (bit_t dn : basis.my_dns_for_ups(idx_up)) {
      coeff_t val = term_action(up, dn);
      vec_out(idx) += val * vec_in(idx);
      ++idx;
    }
    ++idx_up;
  }
}

} // namespace hydra::tj_distributed
