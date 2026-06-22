// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <utility>

#include <xdiag/bits/bitmask.hpp>
#include <xdiag/bits/get_set.hpp>
#include <xdiag/bits/popcount.hpp>
#include <xdiag/bits/zero_one.hpp>
#include <xdiag/kernels/terms/term_offdiag_fermionic.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/logger.hpp>

namespace xdiag::kernels::fermion {

template <typename coeff_t, class basis_t, class fill_f>
void term_hop(Coeff const &c, Op const &op, basis_t const &basis_in,
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

  auto non_zero_term = [&](bit_t const &spins) -> bool {
    return bits::popcount(spins & flipmask) & 1;
  };

  auto term_action = [&](bit_t spins) -> std::pair<bit_t, coeff_t> {
    bool fermi = bits::popcount(spins & fermimask) & 1;
    spins ^= flipmask;
    return {spins, fermi ? t : -t};
  };
  term_offdiag_fermionic(basis_in, basis_out, non_zero_term, term_action, fill);
}
XDIAG_CATCH

} // namespace xdiag::kernels::fermion
