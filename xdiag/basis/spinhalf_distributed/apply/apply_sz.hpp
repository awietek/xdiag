// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <algorithm>
#include <tuple>

#include <xdiag/bits/bitops.hpp>
#include <xdiag/extern/armadillo/armadillo>
#include <xdiag/operators/op.hpp>

namespace xdiag::basis::spinhalf_distributed {

template <class basis_t, typename coeff_t>
void apply_sz(Coupling const &cpl, Op const &op, basis_t const &basis,
              arma::Col<coeff_t> const &vec_in,
              arma::Col<coeff_t> &vec_out) try {
  using bit_t = typename basis_t::bit_t;

  coeff_t H = cpl.scalar().as<coeff_t>();
  int64_t s = op[0];

  bit_t mask = ((bit_t)1 << s);
  coeff_t val_up = H / 2.;
  coeff_t val_dn = -H / 2.;
  int n_postfix_bits = basis.n_postfix_bits();

  int64_t idx = 0;
  for (auto prefix : basis.prefixes()) {
    auto postfixes = basis.postfix_states(prefix);

    // site in postfixes
    if (s < n_postfix_bits) {
      for (auto postfix : postfixes) {
        coeff_t val = ((postfix & mask) ? val_up : val_dn);
        vec_out(idx) += val * vec_in(idx);
        ++idx;
      }

    }
    // site in prefixes
    else {
      bit_t prefix_shifted = (prefix << n_postfix_bits);
      coeff_t val = (prefix_shifted & mask) ? val_up : val_dn;

      for (int64_t start = idx; idx < start + (int64_t)postfixes.size();
           ++idx) {
        vec_out(idx) += val * vec_in(idx);
      }
    }
  } // for (auto prefix : basis.prefixes())
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

} // namespace xdiag::basis::spinhalf_distributed
