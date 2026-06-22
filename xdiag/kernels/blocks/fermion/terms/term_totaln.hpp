// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/bits/popcount.hpp>
#include <xdiag/kernels/terms/term_diag.hpp>
#include <xdiag/utils/error.hpp>

namespace xdiag::kernels::fermion {

// Total particle number operator N_tot = sum_i n_i. It is diagonal in the
// occupation basis with eigenvalue equal to the number of occupied sites,
// which is simply the popcount of the bit pattern.
template <typename coeff_t, class basis_t, class fill_f>
void term_totaln(Coeff const &c, basis_t const &basis_in,
                 basis_t const &basis_out, fill_f fill) try {
  using bit_t = typename basis_t::bit_t;

  coeff_t mu = c.scalar().as<coeff_t>();

  term_diag(
      basis_in, basis_out,
      [&](bit_t spins) -> coeff_t {
        return mu * (coeff_t)bits::popcount(spins);
      },
      fill);
}
XDIAG_CATCH

} // namespace xdiag::kernels::fermion
