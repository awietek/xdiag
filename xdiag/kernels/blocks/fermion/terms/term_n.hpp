// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/bits/zero_one.hpp>
#include <xdiag/bits/get_set.hpp>
#include <xdiag/kernels/terms/term_diag.hpp>
#include <xdiag/utils/error.hpp>

namespace xdiag::kernels::fermion {

template <typename coeff_t, class basis_t, class fill_f>
void term_n(Coeff const &c, Op const &op, basis_t const &basis_in,
            basis_t const &basis_out, fill_f fill) try {
  using bit_t = typename basis_t::bit_t;

  int64_t nsites = basis_in.nsites();
  coeff_t mu = c.scalar().as<coeff_t>();
  int64_t s = op[0];
  bit_t mask = bits::zero<bit_t>(nsites);
  bits::set(mask, s);

  term_diag(
      basis_in, basis_out,
      [&](bit_t spins) -> coeff_t {
        return bits::nonzero(spins & mask) ? mu : 0;
      },
      fill);
}
XDIAG_CATCH

} // namespace xdiag::kernels::fermion
