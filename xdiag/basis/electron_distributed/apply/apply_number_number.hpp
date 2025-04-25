// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/basis/electron_distributed/apply/generic_term_diag.hpp>
#include <xdiag/bits/gbit.hpp>

namespace xdiag::basis::electron_distributed {

template <typename coeff_t, class basis_t>
void apply_ntot_ntot(Coupling const &cpl, Op const &op, basis_t const &basis,
                     const coeff_t *vec_in, coeff_t *vec_out) try {
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

  electron_distributed::generic_term_diag<coeff_t>(basis, apply, vec_in,
                                                   vec_out);

} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <typename coeff_t, class basis_t>
void apply_nupdn(Coupling const &cpl, Op const &op, basis_t const &basis,
                 const coeff_t *vec_in, coeff_t *vec_out) try {
  using bit_t = typename basis_t::bit_t;
  using bits::gbit;

  coeff_t mu = cpl.scalar().as<coeff_t>();
  int64_t s = op[0];
  bit_t mask = (bit_t)1 << s;
  auto apply = [&](bit_t ups, bit_t dns) {
    return (ups & mask & dns) ? mu : 0.;
  };

  electron_distributed::generic_term_diag<coeff_t>(basis, apply, vec_in,
                                                   vec_out);

} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <typename coeff_t, class basis_t>
void apply_nupdn_nupdn(Coupling const &cpl, Op const &op, basis_t const &basis,
                       const coeff_t *vec_in, coeff_t *vec_out) try {
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

  electron_distributed::generic_term_diag<coeff_t>(basis, apply, vec_in,
                                                   vec_out);

} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <typename coeff_t, class basis_t>
void apply_nup_ndn(Coupling const &cpl, Op const &op, basis_t const &basis,
                 fill_f fill) try {
    using bit_t = typename basis_t::bit_t;
    using bits::gbit;

    coeff_t mu = cpl.scalar().as<coeff_t>();
    int64_t s1 = op[0];
    int64_t s2 = op[1];

    bit_t mask1 = (bit_t)1 << s1;
    bit_t mask2 = (bit_t)1 << s2;

    auto apply = [&](bit_t ups, bit_t dns) {
        return (ups & mask1) && (dns & mask2) ? mu : 0.;
    };
    electron_distributed::generic_term_diag<bit_t, coeff_t, symmetric>(cpl, basis, apply, fill);

} catch (Error const &e) {
    XDIAG_RETHROW(e);
}

template <typename coeff_t, class basis_t>
void apply_nup_nup(Coupling const &cpl, Op const &op, basis_t const &basis,
                 fill_f fill) try {
    using bit_t = typename basis_t::bit_t;
    using bits::gbit;

    coeff_t mu = cpl.scalar().as<coeff_t>();
    int64_t s1 = op[0];
    int64_t s2 = op[1];

    bit_t mask1 = (bit_t)1 << s1;
    bit_t mask2 = (bit_t)1 << s2;

    auto apply = [&](bit_t ups, bit_t dns) {
        (void) dns;
        return (ups & mask1) && (ups & mask2) ? 0. : mu;
    };
    electron_distributed::generic_term_diag<bit_t, coeff_t, symmetric>(cpl, basis, apply, fill);

} catch (Error const &e) {
    XDIAG_RETHROW(e);
}

template <typename coeff_t, class basis_t>
void apply_ndn_ndn(Coupling const &cpl, Op const &op, basis_t const &basis,
                 fill_f fill) try {
    using bit_t = typename basis_t::bit_t;
    using bits::gbit;

    coeff_t mu = cpl.scalar().as<coeff_t>();
    int64_t s1 = op[0];
    int64_t s2 = op[1];

    bit_t mask1 = (bit_t)1 << s1;
    bit_t mask2 = (bit_t)1 << s2;

    auto apply = [&](bit_t ups, bit_t dns) {
        (void) ups;
        return (dns & mask1) && (dns & mask2) ? 0. : mu;
    };
    electron_distributed::generic_term_diag<bit_t, coeff_t, symmetric>(cpl, basis, apply, fill);

} catch (Error const &e) {
    XDIAG_RETHROW(e);
}
} // namespace xdiag::basis::electron_distributed
