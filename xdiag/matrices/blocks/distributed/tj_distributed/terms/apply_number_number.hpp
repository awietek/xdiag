// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/matrices/blocks/distributed/tj_distributed/terms/generic_term_diag.hpp>
#include <xdiag/bits/get_set.hpp>

namespace xdiag::basis::tj_distributed {

template <typename coeff_t, class basis_t>
void apply_ntot_ntot(Coeff const &cpl, Op const &op, basis_t const &basis,
                     const coeff_t *vec_in, coeff_t *vec_out) {
  using bit_t = typename basis_t::bit_t;
  using bits::get;

  coeff_t mu = cpl.scalar().as<coeff_t>();
  int64_t s1 = op[0];
  int64_t s2 = op[1];
  auto apply = [&](bit_t ups, bit_t dns) {
    int n1 = get(ups, s1) + get(dns, s1);
    int n2 = get(ups, s2) + get(dns, s2);
    return mu * (coeff_t)(n1 * n2);
  };

  tj_distributed::generic_term_diag<coeff_t>(basis, apply, vec_in, vec_out);
}

} // namespace xdiag::basis::tj_distributed
