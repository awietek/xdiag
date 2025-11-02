// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/basis/spinhalf/apply/apply_term_diag.hpp>
#include <xdiag/bits/bitops.hpp>
#include <xdiag/common.hpp>

namespace xdiag::basis::spinhalf {

// Sz term: H S^z_i
template <typename coeff_t, bool symmetric, class basis_t, class fill_f>
void apply_sz(Coupling const &cpl, Op const &op, basis_t const &basis_in,
              basis_t const &basis_out, fill_f fill) {
  using bit_t = typename basis_t::bit_t;

  coeff_t H = cpl.scalar().as<coeff_t>();
  int64_t s = op[0];
  bit_t mask = ((bit_t)1 << s);

  coeff_t val_up = H / 2.;
  coeff_t val_dn = -H / 2.;

  if (basis_in == basis_out) {

    auto term_coeff = [&mask, &val_up, &val_dn](bit_t spins) -> coeff_t {
      if (spins & mask) {
        return val_up;
      } else {
        return val_dn;
      }
    };

    spinhalf::apply_term_diag<coeff_t>(basis_in, term_coeff, fill);
  } else {
    auto non_zero_term = [](bit_t spins) -> bool {
      return true;
      (void)spins;
    };
    auto term_action = [&](bit_t spins) -> std::pair<bit_t, coeff_t> {
      if (spins & mask) {
        return {spins, val_up};
      } else {
        return {spins, val_dn};
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
