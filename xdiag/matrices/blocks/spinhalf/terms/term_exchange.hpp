// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <utility>

#include <xdiag/bits/get_set.hpp>
#include <xdiag/bits/popcount.hpp>
#include <xdiag/matrices/terms/term_offdiag.hpp>

namespace xdiag::matrices::spinhalf {

template <typename coeff_t, class basis_t, class fill_f>
void term_exchange(Coeff const &c, Op const &op, basis_t const &basis_in,
                   basis_t const &basis_out, fill_f fill) try {
  using bit_t = typename basis_t::bit_t;

  coeff_t J = c.scalar().as<coeff_t>();
  int64_t s1 = op[0];
  int64_t s2 = op[1];

  if (s1 == s2) {
    term_diag(basis_in, basis_out, [&](bit_t) { return J / 2.0; }, fill);
    return;
  }

  bit_t mask = bit_t();
  bits::set(mask, s1);
  bits::set(mask, s2);

  auto non_zero_term = [&](bit_t spins) -> bool {
    return bits::popcount(spins & mask) & 1;
  };

  coeff_t Jhalf = J / 2.0;
  term_offdiag(
      basis_in, basis_out, non_zero_term,
      [&](bit_t spins) -> std::pair<bit_t, coeff_t> {
        return {spins ^ mask, Jhalf};
      },
      fill);
}
XDIAG_CATCH

} // namespace xdiag::matrices::spinhalf
