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

// Nup{s} = n^up_s
template <typename coeff_t, class basis_t, class fill_f>
void term_nup(Coeff const &c, Op const &op,
              basis_t const &basis_in,
              basis_t const &basis_out,
              fill_f fill) try {
  using bit_t = typename basis_t::bit_t;
  coeff_t mu = c.scalar().as<coeff_t>();
  int64_t s = op[0];
  term_diag<coeff_t>(
      basis_in, basis_out,
      [&](bit_t ups, bit_t) -> coeff_t {
        return bits::get(ups, s) ? mu : coeff_t(0);
      },
      fill);
}
XDIAG_CATCH

// Ndn{s} = n^dn_s
template <typename coeff_t, class basis_t, class fill_f>
void term_ndn(Coeff const &c, Op const &op,
              basis_t const &basis_in,
              basis_t const &basis_out,
              fill_f fill) try {
  using bit_t = typename basis_t::bit_t;
  coeff_t mu = c.scalar().as<coeff_t>();
  int64_t s = op[0];
  term_diag<coeff_t>(
      basis_in, basis_out,
      [&](bit_t, bit_t dns) -> coeff_t {
        return bits::get(dns, s) ? mu : coeff_t(0);
      },
      fill);
}
XDIAG_CATCH

// Nupdn{s} = n^up_s n^dn_s (local double occupancy)
template <typename coeff_t, class basis_t, class fill_f>
void term_nupdn(Coeff const &c, Op const &op,
                basis_t const &basis_in,
                basis_t const &basis_out,
                fill_f fill) try {
  using bit_t = typename basis_t::bit_t;
  coeff_t u = c.scalar().as<coeff_t>();
  int64_t s = op[0];
  term_diag<coeff_t>(
      basis_in, basis_out,
      [&](bit_t ups, bit_t dns) -> coeff_t {
        return (bits::get(ups, s) && bits::get(dns, s)) ? u : coeff_t(0);
      },
      fill);
}
XDIAG_CATCH

} // namespace xdiag::matrices::electron
