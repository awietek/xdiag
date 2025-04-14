// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/basis/electron/apply/generic_term_diag.hpp>
#include <xdiag/bits/gbit.hpp>

namespace xdiag::basis::electron {

template <typename coeff_t, bool symmetric, class basis_t, class fill_f>
void apply_ntot_ntot(Coupling const &cpl, Op const &op, basis_t const &basis,
                     fill_f fill) try {
  using bit_t = typename basis_t::bit_t;
  using bits::gbit;

  coeff_t mu = cpl.scalar().as<coeff_t>();
  int64_t s1 = op[0];
  int64_t s2 = op[1];
  auto apply = [&](bit_t ups, bit_t dns) {
    int n1 = gbit(ups, s1) + gbit(dns, s1);
    int n2 = gbit(ups, s2) + gbit(dns, s2);
    return mu * (coeff_t)(n1 * n2);
  };

  generic_term_diag<bit_t, coeff_t, symmetric>(cpl, basis, apply, fill);

} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <typename coeff_t, bool symmetric, class basis_t, class fill_f>
void apply_nupdn(Coupling const &cpl, Op const &op, basis_t const &basis,
                 fill_f fill) try {
  using bit_t = typename basis_t::bit_t;
  using bits::gbit;

  coeff_t mu = cpl.scalar().as<coeff_t>();
  int64_t s = op[0];
  bit_t mask = (bit_t)1 << s;
  auto apply = [&](bit_t ups, bit_t dns) {
    return (ups & mask & dns) ? mu : 0.;
  };

  generic_term_diag<bit_t, coeff_t, symmetric>(cpl, basis, apply, fill);

} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <typename coeff_t, bool symmetric, class basis_t, class fill_f>
void apply_nupdn_nupdn(Coupling const &cpl, Op const &op, basis_t const &basis,
                       fill_f fill) try {
  using bit_t = typename basis_t::bit_t;
  using bits::gbit;

  coeff_t mu = cpl.scalar().as<coeff_t>();
  int64_t s1 = op[0];
  int64_t s2 = op[1];
  bit_t mask1 = (bit_t)1 << s1;
  bit_t mask2 = (bit_t)1 << s2;

  auto apply = [&](bit_t ups, bit_t dns) {
    return (ups & mask1 & dns) && (ups & mask2 & dns) ? mu : 0.;
  };

  generic_term_diag<bit_t, coeff_t, symmetric>(cpl, basis, apply, fill);

} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

} // namespace xdiag::basis::electron
