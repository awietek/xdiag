// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <utility>

#include <xdiag/bits/get_set.hpp>
#include <xdiag/bits/zero_one.hpp>
#include <xdiag/matrices/terms/term_offdiag.hpp>
#include <xdiag/utils/error.hpp>

namespace xdiag::matrices::spinhalf {

// Scalar chirality term: J S*(S x S)

template <typename coeff_t, class basis_t, class fill_f>
void term_scalar_chirality(Coeff const &c, Op const &op,
                           basis_t const &basis_in, basis_t const &basis_out,
                           fill_f fill) try {
  using bit_t = typename basis_t::bit_t;

  complex J = c.scalar().as<complex>();
  coeff_t Jquarter = 0.;
  coeff_t Jquarter_acyclic = 0.;
  if constexpr (isreal<coeff_t>()) {
    XDIAG_THROW("Scalar chirality term cannot be used with real coefficients");
  } else {
    // This realizes ScalarChirality = S_i . (S_j x S_k). The acyclic (inverse)
    // rotation carries -Jquarter, i.e. the operator is (i/4) J (cyclic -
    // acyclic). Using -Jquarter (rather than conj(Jquarter)) keeps the term
    // linear in J for complex coefficients, while still reducing to the
    // Hermitian operator for real J.
    Jquarter = complex(0, 0.25) * J;
    Jquarter_acyclic = -Jquarter;
  }

  bit_t mask = bit_t();
  int64_t s1 = op[0];
  int64_t s2 = op[1];
  int64_t s3 = op[2];
  bits::set(mask, s1);
  bits::set(mask, s2);
  bits::set(mask, s3);

  // scalar chirality annihilates 000 and 111
  auto non_zero_term = [&](bit_t spins) -> bool {
    bit_t threespins = spins & mask;
    return bits::nonzero(threespins) && (threespins != mask);
  };

  // rotate three sites cyclic
  auto term_action_cyclic = [&](bit_t spins) -> std::pair<bit_t, coeff_t> {
    bit_t threespins_cyclic = bit_t();
    bits::set(threespins_cyclic, s2, bits::get(spins, s1));
    bits::set(threespins_cyclic, s3, bits::get(spins, s2));
    bits::set(threespins_cyclic, s1, bits::get(spins, s3));
    bit_t spins_void = spins & (~mask);
    bit_t spins_cyclic = spins_void | threespins_cyclic;
    return {spins_cyclic, Jquarter};
  };

  // rotate three sites acyclic
  auto term_action_acyclic = [&](bit_t spins) -> std::pair<bit_t, coeff_t> {
    bit_t threespins_acyclic = bit_t();
    bits::set(threespins_acyclic, s3, bits::get(spins, s1));
    bits::set(threespins_acyclic, s1, bits::get(spins, s2));
    bits::set(threespins_acyclic, s2, bits::get(spins, s3));
    bit_t spins_void = spins & (~mask);
    bit_t spins_acyclic = spins_void | threespins_acyclic;
    return {spins_acyclic, Jquarter_acyclic};
  };
  term_offdiag(basis_in, basis_out, non_zero_term, term_action_cyclic, fill);
  term_offdiag(basis_in, basis_out, non_zero_term, term_action_acyclic, fill);
}
XDIAG_CATCH

} // namespace xdiag::matrices::spinhalf
