// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/basis/tj/apply/generic_term_diag.hpp>
#include <xdiag/bits/gbit.hpp>
#include <xdiag/common.hpp>
#include <xdiag/operators/op.hpp>

namespace xdiag::basis::tj {

template <bool symmetric, typename coeff_t, typename basis_t, typename fill_f>
void apply_number_number(Coupling const &cpl, Op const &op,
                         basis_t const &basis_in, basis_t const &basis_out,
                         fill_f fill) {
  using bit_t = typename basis_t::bit_t;
  using bits::gbit;
  coeff_t mu = cpl.scalar().as<coeff_t>();
  int64_t s1 = op[0];
  int64_t s2 = op[1];

  auto term_action = [&](bit_t up, bit_t dn) {
    int n1 = gbit(up, s1) + gbit(dn, s1);
    int n2 = gbit(up, s2) + gbit(dn, s2);
    return mu * (coeff_t)(n1 * n2);
  };
  tj::generic_term_diag<symmetric, coeff_t>(basis_in, basis_out, term_action,
                                            fill);
}

} // namespace xdiag::basis::tj
