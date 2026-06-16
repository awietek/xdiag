// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <utility>

#include <xdiag/bits/get_set.hpp>
#include <xdiag/bits/zero_one.hpp>
#include <xdiag/matrices/terms/term_offdiag.hpp>
#include <xdiag/utils/error.hpp>

namespace xdiag::matrices::spinhalf {

template <typename coeff_t, class basis_t, class fill_f>
void term_spsm(Coeff const &c, Op const &op, basis_t const &basis_in,
               basis_t const &basis_out, fill_f fill) try {
  using bit_t = typename basis_t::bit_t;

  coeff_t J = c.scalar().as<coeff_t>();
  int64_t s = op[0];
  bit_t mask = bit_t();
  bits::set(mask, s);

  if (op.type() == "S+") {
    auto non_zero_term = [&](bit_t spins) {
      return !bits::nonzero(spins & mask);
    };
    auto term_action = [&](bit_t spins) -> std::pair<bit_t, coeff_t> {
      bit_t spins_flip = spins | mask;
      return {spins_flip, J};
    };
    term_offdiag(basis_in, basis_out, non_zero_term, term_action, fill);
  } else { // op.type() == "S-"
    auto non_zero_term = [&](bit_t spins) {
      return bits::nonzero(spins & mask);
    };
    auto term_action = [&](bit_t spins) -> std::pair<bit_t, coeff_t> {
      bit_t spins_flip = spins ^ mask;
      return {spins_flip, J};
    };
    term_offdiag(basis_in, basis_out, non_zero_term, term_action, fill);
  }
}
XDIAG_CATCH

} // namespace xdiag::matrices::spinhalf
