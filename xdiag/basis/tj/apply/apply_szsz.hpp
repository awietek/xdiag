// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/basis/tj/apply/generic_term_diag.hpp>
#include <xdiag/common.hpp>
#include <xdiag/operators/op.hpp>

namespace xdiag::basis::tj {

template <typename coeff_t, bool symmetric, class basis_t, class fill_f>
void apply_szsz(Coupling const &cpl, Op const &op, basis_t const &basis,
                fill_f fill) try {
  using bit_t = typename basis_t::bit_t;

  coeff_t J = cpl.scalar().as<coeff_t>();
  int64_t s1 = op[0];
  int64_t s2 = op[1];
  bit_t s1_mask = (bit_t)1 << s1;
  bit_t s2_mask = (bit_t)1 << s2;

  // Set values for same/diff (tJ block definition)
  coeff_t val_same, val_diff;
  if (op.type() == "SzSz") {
    val_same = J / 4.;
    val_diff = -J / 4.;
  } else { // (type == "tJSzSz")
    val_same = 0.;
    val_diff = -J / 2.;
  }

  auto term_action = [&](bit_t up, bit_t dn) -> coeff_t {
    bool b1_up = up & s1_mask;
    bool b2_up = up & s2_mask;
    bool b1_dn = dn & s1_mask;
    bool b2_dn = dn & s2_mask;

    if ((b1_up && b2_up) || (b1_dn && b2_dn)) {
      return val_same;
    } else if ((b1_up && b2_dn) || (b1_dn && b2_up)) {
      return val_diff;
    } else {
      return 0.;
    }
  };

  tj::generic_term_diag<coeff_t, symmetric>(basis, term_action, fill);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

} // namespace xdiag::basis::tj
