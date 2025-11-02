// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/basis/spinhalf/apply/apply_term_diag.hpp>
#include <xdiag/bits/bitops.hpp>
#include <xdiag/common.hpp>

namespace xdiag::basis::spinhalf {

// Ising term: J S^z_i S^z_j

template <typename coeff_t, bool symmetric, class basis_t, class fill_f>
void apply_szsz(Coupling const &cpl, Op const &op, basis_t const &basis_in,
                basis_t const &basis_out, fill_f fill) {
  using bit_t = typename basis_t::bit_t;

  coeff_t J = cpl.scalar().as<coeff_t>();
  int64_t s1 = op[0];
  int64_t s2 = op[1];
  bit_t mask = ((bit_t)1 << s1) | ((bit_t)1 << s2);

  coeff_t val_same = J / 4.0;
  coeff_t val_diff = -J / 4.0;

  if (basis_in == basis_out) {

    auto term_coeff = [&](bit_t spins) -> coeff_t {
      if (bits::popcnt(spins & mask) & 1) {
        return val_diff;
      } else {
        return val_same;
      }
    };
    spinhalf::apply_term_diag<coeff_t>(basis_in, term_coeff, fill);

  } else {
    auto non_zero_term = [](bit_t spins) -> bool {
      return true;
      (void)spins;
    };
    auto term_action = [&](bit_t spins) -> std::pair<bit_t, coeff_t> {
      if (bits::popcnt(spins & mask) & 1) {
        return {spins, val_diff};
      } else {
        return {spins, val_same};
      }
    };
    if constexpr (symmetric) {
      spinhalf::apply_term_offdiag_sym<coeff_t>(
          basis_in, basis_out, non_zero_term, term_action, fill);
    } else {
      spinhalf::apply_term_offdiag_no_sym<coeff_t>(
          basis_in, basis_out, non_zero_term, term_action, fill);
    }
  }
}

} // namespace xdiag::basis::spinhalf
