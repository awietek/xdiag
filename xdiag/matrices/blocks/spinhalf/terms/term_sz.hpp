// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/bits/get_set.hpp>
#include <xdiag/bits/zero_one.hpp>
#include <xdiag/matrices/terms/term_diag.hpp>
#include <xdiag/utils/error.hpp>

namespace xdiag::matrices::spinhalf {

template <typename coeff_t, class basis_t, class fill_f>
void term_sz(Coeff const &c, Op const &op, basis_t const &basis_in,
             basis_t const &basis_out, fill_f fill) try {
  using bit_t = typename basis_t::bit_t;

  coeff_t H = c.scalar().as<coeff_t>();
  int64_t s = op[0];
  bit_t mask = bit_t();
  bits::set(mask, s);

  coeff_t val_up = H / 2.;
  coeff_t val_dn = -H / 2.;

  term_diag(
      basis_in, basis_out,
      [&](bit_t spins) -> coeff_t {
        return bits::nonzero(spins & mask) ? val_up : val_dn;
      },
      fill);
}
XDIAG_CATCH

} // namespace xdiag::matrices::spinhalf
