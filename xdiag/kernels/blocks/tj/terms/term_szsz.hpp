// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>

#include <xdiag/basis/basis_tj.hpp>
#include <xdiag/bits/bitmask.hpp>
#include <xdiag/bits/get_set.hpp>
#include <xdiag/bits/popcount.hpp>
#include <xdiag/kernels/blocks/tj/terms/term_diag.hpp>
#include <xdiag/operators/coeff.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/utils/error.hpp>

namespace xdiag::kernels::tj {

// Diagonal spin-spin coupling on a bond {i,j}. Two conventions (matching the old
// tJ block), selected by the operator type:
//
//   SzSz{i,j}   = Sz_i Sz_j ,  Sz_s = (n^up_s - n^dn_s)/2
//                 -> +J/4 (parallel spins), -J/4 (antiparallel), 0 (a site empty)
//   tJSzSz{i,j} = Sz_i Sz_j - (1/4) n_i n_j  (diagonal part of J (S_i.S_j -
//                 n_i n_j /4), the t-J Heisenberg coupling)
//                 -> 0 (parallel), -J/2 (antiparallel), 0 (a site empty)
//
// Only the value on occupied bonds differs; both vanish if either site is empty.
template <typename coeff_t, class enumeration_t, class fill_f>
void term_szsz(Coeff const &c, Op const &op,
               basis::BasistJ<enumeration_t> const &basis_in,
               basis::BasistJ<enumeration_t> const &basis_out,
               fill_f fill) try {
  using bit_t = typename enumeration_t::bit_t;
  coeff_t J = c.scalar().as<coeff_t>();
  int64_t i = op[0];
  int64_t j = op[1];
  int64_t nsites = basis_in.nsites();
  bit_t sitesmask = bits::bitmask<bit_t>(nsites, nsites);

  coeff_t val_same, val_diff;
  if (op.type() == "tJSzSz") {
    val_same = coeff_t(0.0);
    val_diff = -J / coeff_t(2.0);
  } else { // "SzSz"
    val_same = J / coeff_t(4.0);
    val_diff = -J / coeff_t(4.0);
  }

  term_diag<coeff_t>(
      basis_in, basis_out,
      [=](bit_t ups) {
        bool up_i = bits::get(ups, i);
        bool up_j = bits::get(ups, j);
        int64_t ri = tj_dn_rank<bit_t>(ups, i, nsites, sitesmask);
        int64_t rj = tj_dn_rank<bit_t>(ups, j, nsites, sitesmask);
        return [=](bit_t dnc) -> coeff_t {
          bool dn_i = (ri >= 0) && bits::get(dnc, ri);
          bool dn_j = (rj >= 0) && bits::get(dnc, rj);
          if ((up_i && up_j) || (dn_i && dn_j)) {
            return val_same;
          } else if ((up_i && dn_j) || (dn_i && up_j)) {
            return val_diff;
          } else {
            return coeff_t(0.0);
          }
        };
      },
      fill);
}
XDIAG_CATCH

} // namespace xdiag::kernels::tj
