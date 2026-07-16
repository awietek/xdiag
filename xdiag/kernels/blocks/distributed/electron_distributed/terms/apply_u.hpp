// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/kernels/blocks/distributed/electron_distributed/terms/generic_term_diag.hpp>
#include <xdiag/math/complex.hpp>
#include <xdiag/operators/coeff.hpp>
#include <xdiag/operators/op.hpp>

namespace xdiag::basis::electron_distributed {

template <typename coeff_t, class basis_t>
void apply_u(Coeff const &cpl, basis_t const &basis, const coeff_t *vec_in,
             coeff_t *vec_out) {
  using bit_t = typename basis_t::bit_t;

  coeff_t U = cpl.scalar().as<coeff_t>();
  auto term_action = [&](bit_t up, bit_t dn) {
    return U * (double)bits::popcount(up & dn);
  };
  electron_distributed::generic_term_diag<coeff_t>(basis, term_action, vec_in,
                                                   vec_out);
}

} // namespace xdiag::basis::electron_distributed
