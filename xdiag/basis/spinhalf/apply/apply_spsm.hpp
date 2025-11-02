// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <functional>
#include <string>

#include <xdiag/basis/spinhalf/apply/apply_term_offdiag_no_sym.hpp>
#include <xdiag/basis/spinhalf/apply/apply_term_offdiag_sym.hpp>
#include <xdiag/bits/bitops.hpp>
#include <xdiag/common.hpp>

namespace xdiag::basis::spinhalf {

// S+ or S- term: J S^+_i   OR   J S^-_i

template <typename coeff_t, bool symmetric, class basis_t, class fill_f>
void apply_spsm(Coupling const &cpl, Op const &op, basis_t const &basis_in,
                basis_t const &basis_out, fill_f fill) {
  using bit_t = typename basis_t::bit_t;

  coeff_t J = cpl.scalar().as<coeff_t>();
  int64_t s = op[0];
  bit_t mask = ((bit_t)1 << s);

  // Define actions of op
  std::function<bool(bit_t)> non_zero_term;
  std::function<std::pair<bit_t, coeff_t>(bit_t)> term_action;

  if (op.type() == "S+") {
    non_zero_term = [&](bit_t spins) { return !(spins & mask); };
    term_action = [&](bit_t spins) -> std::pair<bit_t, coeff_t> {
      bit_t spins_flip = spins | mask;
      return {spins_flip, J};
    };
  } else { // op.type() == "S-"
    non_zero_term = [&](bit_t spins) { return spins & mask; };
    term_action = [&](bit_t spins) -> std::pair<bit_t, coeff_t> {
      bit_t spins_flip = spins ^ mask;
      return {spins_flip, J};
    };
  }

  // Dispatch either symmetric of unsymmetric term application
  if constexpr (symmetric) {
    spinhalf::apply_term_offdiag_sym<coeff_t>(basis_in, basis_out,
                                              non_zero_term, term_action, fill);
  } else {
    spinhalf::apply_term_offdiag_no_sym<coeff_t>(
        basis_in, basis_out, non_zero_term, term_action, fill);
  }
}

} // namespace xdiag::basis::spinhalf
