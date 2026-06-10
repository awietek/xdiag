// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <utility>

#include <xdiag/bits/bitmask.hpp>
#include <xdiag/bits/get_set.hpp>
#include <xdiag/bits/nonzero.hpp>
#include <xdiag/bits/popcount.hpp>
#include <xdiag/bits/zero_one.hpp>
#include <xdiag/matrices/terms/term_offdiag.hpp>
#include <xdiag/utils/error.hpp>

namespace xdiag::matrices::fermion {

template <typename coeff_t, class basis_t, class fill_f>
void term_hop_asym(Coeff const &c, Op const &op, basis_t const &basis_in,
                   basis_t const &basis_out, fill_f fill) try {
  using bit_t = typename basis_t::bit_t;

  int64_t nsites = basis_in.nsites();
  coeff_t t = c.scalar().as<coeff_t>();
  int64_t s1 = op[0];
  int64_t s2 = op[1];
  bit_t flipmask = bits::zero<bit_t>(nsites);
  bits::set(flipmask, s1);
  bits::set(flipmask, s2);
  int64_t l = std::min(s1, s2);
  int64_t u = std::max(s1, s2);
  bit_t fermimask = bits::bitmask<bit_t>(nsites, u - l - 1) << (l + 1);

  auto non_zero_term = [&flipmask](bit_t const &spins) -> bool {
    return bits::popcount(spins & flipmask) & 1;
  };

  auto term_action = [&](bit_t spins) -> std::pair<bit_t, coeff_t> {
    bool fermi = bits::popcount(spins & fermimask) & 1;
    spins ^= flipmask;
    if (bits::get(spins, s2)) {
      return {spins, fermi ? t : -t};
    } else {
      return {spins, fermi ? -t : t};
    }
  };
}
XDIAG_CATCH

} // namespace xdiag::matrices::fermion
