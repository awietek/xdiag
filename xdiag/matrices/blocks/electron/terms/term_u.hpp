// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/bits/popcount.hpp>
#include <xdiag/matrices/blocks/electron/terms/term_diag.hpp>
#include <xdiag/operators/coeff.hpp>
#include <xdiag/utils/error.hpp>

namespace xdiag::matrices::electron {

// HubbardU = U * sum_i n^up_i n^dn_i. Site-free: the number of doubly occupied
// sites of a product state (ups, dns) is popcount(ups & dns), so the diagonal
// entry is U times that count.
template <typename coeff_t, class basis_t, class fill_f>
void term_hubbardu(Coeff const &c,
                   basis_t const &basis_in,
                   basis_t const &basis_out,
                   fill_f fill) try {
  using bit_t = typename basis_t::bit_t;
  coeff_t U = c.scalar().as<coeff_t>();
  term_diag<coeff_t>(
      basis_in, basis_out,
      [&](bit_t ups, bit_t dns) -> coeff_t {
        return U * (coeff_t)bits::popcount(ups & dns);
      },
      fill);
}
XDIAG_CATCH

} // namespace xdiag::matrices::electron
