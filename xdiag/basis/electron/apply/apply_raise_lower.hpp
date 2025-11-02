// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <string>

#include <xdiag/basis/electron/apply/generic_term_dns.hpp>
#include <xdiag/basis/electron/apply/generic_term_ups.hpp>
#include <xdiag/bits/bitops.hpp>
#include <xdiag/parallel/omp/omp_utils.hpp>

namespace xdiag::basis::electron {

template <bool symmetric, typename coeff_t, typename basis_t, typename fill_f>
void apply_raise_lower(Coupling const &cpl, Op const &op,
                       basis_t const &basis_in, basis_t const &basis_out,
                       fill_f fill) {
  using bit_t = typename basis_t::bit_t;

  coeff_t c = cpl.scalar().as<coeff_t>();
  int64_t s = op[0];
  bit_t site_mask = (bit_t)1 << s;
  bit_t fermi_mask = site_mask - 1;
  std::string type = op.type();

  // Raising operators
  if ((type == "Cdagup") || (type == "Cdagdn")) {
    auto non_zero_term = [&](bit_t const &spins) -> bool {
      return (spins & site_mask) == 0;
    };
    auto term_action = [&](bit_t spins) -> std::pair<bit_t, coeff_t> {
      bool fermi = bits::popcnt(spins & fermi_mask) & 1;
      return {spins ^ site_mask, fermi ? -c : c};
    };

    if (type == "Cdagup") {
      electron::generic_term_ups<symmetric, coeff_t>(
          basis_in, basis_out, non_zero_term, term_action, fill);
    } else if (type == "Cdagdn") {
      electron::generic_term_dns<symmetric, coeff_t, true>(
          basis_in, basis_out, non_zero_term, term_action, fill);
    }

    // Lowering operators
  } else if ((type == "Cup") || (type == "Cdn")) {
    auto non_zero_term = [&](bit_t const &spins) -> bool {
      return (spins & site_mask);
    };
    auto term_action = [&](bit_t spins) -> std::pair<bit_t, coeff_t> {
      bool fermi = bits::popcnt(spins & fermi_mask) & 1;
      return {spins ^ site_mask, fermi ? -c : c};
    };

    if (type == "Cup") {
      electron::generic_term_ups<symmetric, coeff_t>(
          basis_in, basis_out, non_zero_term, term_action, fill);
    } else if (type == "Cdn") {
      electron::generic_term_dns<symmetric, coeff_t, true>(
          basis_in, basis_out, non_zero_term, term_action, fill);
    }
  }
}
} // namespace xdiag::basis::electron
