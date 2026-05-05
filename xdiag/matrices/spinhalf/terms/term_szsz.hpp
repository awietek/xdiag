// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/bits/get_set_bit.hpp>
#include <xdiag/bits/popcount.hpp>
#include <xdiag/matrices/terms/term_diag.hpp>
#include <xdiag/utils/error.hpp>

namespace xdiag::matrices::spinhalf {

template <typename coeff_t, class basis_t, class fill_f>
void term_szsz(Coeff const &c, Op const &op, basis_t const &basis_in,
               basis_t const &basis_out, fill_f fill) try {
  using bit_t = typename basis_t::bit_t;

  coeff_t J = c.scalar().as<coeff_t>();

  if (op[0] == op[1]) {
    term_diag(basis_in, basis_out, [&](bit_t) { return J / 4.0; }, fill);
    return;
  }

  coeff_t val_same = J / 4.0;
  coeff_t val_diff = -J / 4.0;

  bit_t mask = bit_t();
  bits::set_bit(mask, op[0]);
  bits::set_bit(mask, op[1]);

  term_diag(
      basis_in, basis_out,
      [&](bit_t spins) {
        return bits::popcount(spins & mask) & 1 ? val_diff : val_same;
      },
      fill);
}
XDIAG_CATCH

} // namespace xdiag::matrices::spinhalf
