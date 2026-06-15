// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <utility>

#include <xdiag/matrices/terms/cdagc_string.hpp>
#include <xdiag/matrices/terms/term_offdiag_fermionic.hpp>
#include <xdiag/operators/coeff.hpp>
#include <xdiag/operators/monomial.hpp>
#include <xdiag/utils/error.hpp>

namespace xdiag::matrices::fermion {

// Applies a normal-ordered string of elementary fermion operators (Cdag/C) to a
// fermion basis. The string decode and the per-state Jordan-Wigner sign are
// handled by the shared CdagCString (see matrices/terms/cdagc_string.hpp).
template <typename coeff_t, class basis_t, class fill_f>
void term_cdagc_string(Coeff const &c, Monomial const &mono,
                       basis_t const &basis_in, basis_t const &basis_out,
                       fill_f fill) try {
  using bit_t = typename basis_t::bit_t;

  coeff_t cf = c.scalar().as<coeff_t>();
  CdagCString<bit_t> str(basis_in.nsites(), mono, "Cdag", "C");

  auto non_zero_term = [&](bit_t spins) -> bool {
    return str.non_zero(spins);
  };
  auto term_action = [&](bit_t spins) -> std::pair<bit_t, coeff_t> {
    std::pair<bit_t, bool> a = str.action(spins);
    return {a.first, a.second ? -cf : cf};
  };
  term_offdiag_fermionic(basis_in, basis_out, non_zero_term, term_action, fill);
}
XDIAG_CATCH

} // namespace xdiag::matrices::fermion
