// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/basis/spinhalf/apply/apply_term_offdiag_no_sym.hpp>
#include <xdiag/basis/spinhalf/apply/apply_term_offdiag_sym.hpp>
#include <xdiag/bits/bitops.hpp>
#include <xdiag/common.hpp>

namespace xdiag::basis::spinhalf {

// Scalar chirality term: J S*(S x S)

template <typename coeff_t, bool symmetric, class basis_t, class fill_f>
void apply_scalar_chirality(Coupling const &cpl, Op const &op,
                            basis_t const &basis_in, basis_t const &basis_out,
                            fill_f fill) try {
  using bit_t = typename basis_t::bit_t;
  using bits::gbit;

  complex J = cpl.scalar().as<complex>();
  coeff_t Jquarter = 0.;
  coeff_t Jquarter_conj = 0.;
  if constexpr (isreal<coeff_t>()) {
    XDIAG_THROW(
        "Error in spinhalf::apply_scalar_chirality: scalar chirality term "
        "cannot be used with real coefficients");
  } else {
    Jquarter = complex(0, -0.25) * J;
    Jquarter_conj = xdiag::conj(Jquarter);
  }

  int64_t s1 = op[0];
  int64_t s2 = op[1];
  int64_t s3 = op[2];
  bit_t spinmask = ((bit_t)1 << s1) | ((bit_t)1 << s2) | ((bit_t)1 << s3);

  // scalar chirality annihilates 000 and 111
  auto non_zero_term = [&spinmask](bit_t spins) -> bool {
    bit_t threespins = spins & spinmask;
    return (threespins != 0) && (threespins != spinmask);
  };

  // rotate three sites cyclic
  auto term_action_cyclic = [&spinmask, &s1, &s2, &s3, &Jquarter](
                                bit_t spins) -> std::pair<bit_t, coeff_t> {
    bit_t b1 = gbit(spins, s1);
    bit_t b2 = gbit(spins, s2);
    bit_t b3 = gbit(spins, s3);
    bit_t threespins_cyclic = (b1 << s2) | (b2 << s3) | (b3 << s1);
    bit_t spins_void = spins & (~spinmask);
    bit_t spins_cyclic = spins_void | threespins_cyclic;
    return {spins_cyclic, Jquarter};
  };

  // rotate three sites acyclic
  auto term_action_acyclic = [&spinmask, &s1, &s2, &s3, &Jquarter_conj](
                                 bit_t spins) -> std::pair<bit_t, coeff_t> {
    bit_t b1 = gbit(spins, s1);
    bit_t b2 = gbit(spins, s2);
    bit_t b3 = gbit(spins, s3);
    bit_t threespins_acyclic = (b1 << s3) | (b2 << s1) | (b3 << s2);
    bit_t spins_void = spins & (~spinmask);
    bit_t spins_acyclic = spins_void | threespins_acyclic;
    return {spins_acyclic, Jquarter_conj};
  };

  // Dispatch either symmetric of unsymmetric term application
  if constexpr (symmetric) {
    spinhalf::apply_term_offdiag_sym<coeff_t>(
        basis_in, basis_out, non_zero_term, term_action_cyclic, fill);
    spinhalf::apply_term_offdiag_sym<coeff_t>(
        basis_in, basis_out, non_zero_term, term_action_acyclic, fill);

  } else {
    spinhalf::apply_term_offdiag_no_sym<coeff_t>(
        basis_in, basis_out, non_zero_term, term_action_cyclic, fill);
    spinhalf::apply_term_offdiag_no_sym<coeff_t>(
        basis_in, basis_out, non_zero_term, term_action_acyclic, fill);
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

} // namespace xdiag::basis::spinhalf
