// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>

#include <xdiag/bits/get_set.hpp>
#include <xdiag/matrices/blocks/electron/terms/term_diag.hpp>
#include <xdiag/operators/coeff.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/utils/error.hpp>

namespace xdiag::matrices::electron {

// SzSz{i,j} = Sz_i Sz_j with Sz_s = (n^up_s - n^dn_s) / 2.
template <typename coeff_t, class enumeration_t, class fill_f>
void term_szsz(Coeff const &c, Op const &op,
               basis::BasisElectron<enumeration_t> const &basis_in,
               basis::BasisElectron<enumeration_t> const &basis_out,
               fill_f fill) try {
  using bit_t = typename enumeration_t::bit_t;
  coeff_t J = c.scalar().as<coeff_t>();
  int64_t i = op[0], j = op[1];
  term_diag<coeff_t>(
      basis_in, basis_out,
      [&](bit_t ups, bit_t dns) -> coeff_t {
        int64_t two_sz_i = (int64_t)bits::get(ups, i) - (int64_t)bits::get(dns, i);
        int64_t two_sz_j = (int64_t)bits::get(ups, j) - (int64_t)bits::get(dns, j);
        return J * (coeff_t)(two_sz_i * two_sz_j) * coeff_t(0.25);
      },
      fill);
}
XDIAG_CATCH

} // namespace xdiag::matrices::electron
