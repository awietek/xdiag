// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/basis/tj/apply/generic_term_diag.hpp>
#include <xdiag/common.hpp>
#include <xdiag/operators/op.hpp>

namespace xdiag::basis::tj {

template <bool symmetric, typename coeff_t, typename basis_t, typename fill_f>
void apply_number(Coupling const &cpl, Op const &op, basis_t const &basis_in,
                  basis_t const &basis_out, fill_f fill) {
  using bit_t = typename basis_t::bit_t;
  coeff_t mu = cpl.scalar().as<coeff_t>();
  int64_t s = op[0];
  bit_t mask = (bit_t)1 << s;

  if (op.type() == "Nup") {
    auto term_action = [&](bit_t up, bit_t dn) {
      (void)dn;
      return (up & mask) ? mu : 0.;
    };
    tj::generic_term_diag<symmetric, coeff_t>(basis_in, basis_out, term_action,
                                              fill);
  } else { // type == "Ndn"
    auto term_action = [&](bit_t up, bit_t dn) {
      (void)up;
      return (dn & mask) ? mu : 0.;
    };
    tj::generic_term_diag<symmetric, coeff_t>(basis_in, basis_out, term_action,
                                              fill);
  }
}

} // namespace xdiag::basis::tj
