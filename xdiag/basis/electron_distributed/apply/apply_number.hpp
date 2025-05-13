// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/basis/electron_distributed/apply/generic_term_diag.hpp>
#include <xdiag/common.hpp>
#include <xdiag/operators/op.hpp>

namespace xdiag::basis::electron_distributed {

template <typename coeff_t, class basis_t>
void apply_number(Coupling const &cpl, Op const &op, basis_t const &basis,
                  const coeff_t *vec_in, coeff_t *vec_out) {
  using bit_t = typename basis_t::bit_t;

  coeff_t mu = cpl.scalar().as<coeff_t>();
  int64_t s = op[0];
  std::string type = op.type();

  bit_t mask = (bit_t)1 << s;
  if (type == "Nup") {
    auto term_action = [&](bit_t up, bit_t dn) {
      (void)dn;
      return (up & mask) ? mu : 0.;
    };
    electron_distributed::generic_term_diag<coeff_t>(basis, term_action, vec_in,
                                                     vec_out);
  } else if (type == "Ndn") {
    auto term_action = [&](bit_t up, bit_t dn) {
      (void)up;
      return (dn & mask) ? mu : 0.;
    };
    electron_distributed::generic_term_diag<coeff_t>(basis, term_action, vec_in,
                                                     vec_out);
  }
}

} // namespace xdiag::basis::electron_distributed
