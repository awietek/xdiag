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

// NupNup{i,j} = n^up_i n^up_j
template <typename coeff_t, class enumeration_t, class fill_f>
void term_nupnup(Coeff const &c, Op const &op,
                 basis::BasisElectron<enumeration_t> const &basis_in,
                 basis::BasisElectron<enumeration_t> const &basis_out,
                 fill_f fill) try {
  using bit_t = typename enumeration_t::bit_t;
  coeff_t v = c.scalar().as<coeff_t>();
  int64_t i = op[0], j = op[1];
  term_diag<coeff_t>(
      basis_in, basis_out,
      [&](bit_t ups, bit_t) -> coeff_t {
        return (bits::get(ups, i) && bits::get(ups, j)) ? v : coeff_t(0);
      },
      fill);
}
XDIAG_CATCH

// NdnNdn{i,j} = n^dn_i n^dn_j
template <typename coeff_t, class enumeration_t, class fill_f>
void term_ndnndn(Coeff const &c, Op const &op,
                 basis::BasisElectron<enumeration_t> const &basis_in,
                 basis::BasisElectron<enumeration_t> const &basis_out,
                 fill_f fill) try {
  using bit_t = typename enumeration_t::bit_t;
  coeff_t v = c.scalar().as<coeff_t>();
  int64_t i = op[0], j = op[1];
  term_diag<coeff_t>(
      basis_in, basis_out,
      [&](bit_t, bit_t dns) -> coeff_t {
        return (bits::get(dns, i) && bits::get(dns, j)) ? v : coeff_t(0);
      },
      fill);
}
XDIAG_CATCH

// NupNdn{i,j} = n^up_i n^dn_j
template <typename coeff_t, class enumeration_t, class fill_f>
void term_nupndn(Coeff const &c, Op const &op,
                 basis::BasisElectron<enumeration_t> const &basis_in,
                 basis::BasisElectron<enumeration_t> const &basis_out,
                 fill_f fill) try {
  using bit_t = typename enumeration_t::bit_t;
  coeff_t v = c.scalar().as<coeff_t>();
  int64_t i = op[0], j = op[1];
  term_diag<coeff_t>(
      basis_in, basis_out,
      [&](bit_t ups, bit_t dns) -> coeff_t {
        return (bits::get(ups, i) && bits::get(dns, j)) ? v : coeff_t(0);
      },
      fill);
}
XDIAG_CATCH

// NdnNup{i,j} = n^dn_i n^up_j
template <typename coeff_t, class enumeration_t, class fill_f>
void term_ndnnup(Coeff const &c, Op const &op,
                 basis::BasisElectron<enumeration_t> const &basis_in,
                 basis::BasisElectron<enumeration_t> const &basis_out,
                 fill_f fill) try {
  using bit_t = typename enumeration_t::bit_t;
  coeff_t v = c.scalar().as<coeff_t>();
  int64_t i = op[0], j = op[1];
  term_diag<coeff_t>(
      basis_in, basis_out,
      [&](bit_t ups, bit_t dns) -> coeff_t {
        return (bits::get(dns, i) && bits::get(ups, j)) ? v : coeff_t(0);
      },
      fill);
}
XDIAG_CATCH

// NtotNtot{i,j} = (n^up_i + n^dn_i)(n^up_j + n^dn_j)
template <typename coeff_t, class enumeration_t, class fill_f>
void term_ntotntot(Coeff const &c, Op const &op,
                   basis::BasisElectron<enumeration_t> const &basis_in,
                   basis::BasisElectron<enumeration_t> const &basis_out,
                   fill_f fill) try {
  using bit_t = typename enumeration_t::bit_t;
  coeff_t v = c.scalar().as<coeff_t>();
  int64_t i = op[0], j = op[1];
  term_diag<coeff_t>(
      basis_in, basis_out,
      [&](bit_t ups, bit_t dns) -> coeff_t {
        int64_t ntot_i = bits::get(ups, i) + bits::get(dns, i);
        int64_t ntot_j = bits::get(ups, j) + bits::get(dns, j);
        return v * (coeff_t)(ntot_i * ntot_j);
      },
      fill);
}
XDIAG_CATCH

// NupdnNupdn{i,j} = (n^up_i n^dn_i)(n^up_j n^dn_j)
template <typename coeff_t, class enumeration_t, class fill_f>
void term_nupdnnupdn(Coeff const &c, Op const &op,
                     basis::BasisElectron<enumeration_t> const &basis_in,
                     basis::BasisElectron<enumeration_t> const &basis_out,
                     fill_f fill) try {
  using bit_t = typename enumeration_t::bit_t;
  coeff_t v = c.scalar().as<coeff_t>();
  int64_t i = op[0], j = op[1];
  term_diag<coeff_t>(
      basis_in, basis_out,
      [&](bit_t ups, bit_t dns) -> coeff_t {
        bool di = bits::get(ups, i) && bits::get(dns, i);
        bool dj = bits::get(ups, j) && bits::get(dns, j);
        return (di && dj) ? v : coeff_t(0);
      },
      fill);
}
XDIAG_CATCH

} // namespace xdiag::matrices::electron
