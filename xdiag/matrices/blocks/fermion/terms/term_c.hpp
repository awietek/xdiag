// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <utility>

#include <xdiag/bits/zero_one.hpp>
#include <xdiag/bits/bitmask.hpp>
#include <xdiag/bits/get_set.hpp>
#include <xdiag/bits/nonzero.hpp>
#include <xdiag/bits/popcount.hpp>
#include <xdiag/matrices/terms/term_offdiag.hpp>
#include <xdiag/utils/error.hpp>

namespace xdiag::matrices::fermion {

template <typename coeff_t, class basis_t, class fill_f>
void term_c(Coeff const &c, Op const &op, basis_t const &basis_in,
            basis_t const &basis_out, fill_f fill) try {
  using bit_t = typename basis_t::bit_t;

  int64_t nsites = basis_in.nsites();
  coeff_t cf = c.scalar().as<coeff_t>();
  int64_t s = op[0];
  bit_t mask = bits::zero<bit_t>(nsites);
  bits::set(mask, s);
  bit_t fermi_mask = bits::bitmask<bit_t>(nsites, s);

  auto non_zero_term = [&](bit_t spins) { return bits::nonzero(spins & mask); };
  auto term_action = [&](bit_t spins) -> std::pair<bit_t, coeff_t> {
    bit_t spins_flip = spins ^ mask;
    bool fermi = bits::popcount(spins & fermi_mask) & 1;
    return {spins_flip, fermi ? -cf : cf};
  };
  term_offdiag(basis_in, basis_out, non_zero_term, term_action, fill);
}
XDIAG_CATCH

} // namespace xdiag::matrices::fermion
