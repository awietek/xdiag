// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/basis/tj/apply/generic_term_mixed.hpp>
#include <xdiag/bits/bitops.hpp>
#include <xdiag/common.hpp>
#include <xdiag/operators/op.hpp>

namespace xdiag::basis::tj {

template <bool symmetric, typename coeff_t, typename basis_t, typename fill_f>
void apply_exchange(Coupling const &cpl, Op const &op, basis_t const &basis_in,
                    basis_t const &basis_out, fill_f fill) {
  using bit_t = typename basis_t::bit_t;

  coeff_t J = cpl.scalar().as<coeff_t>();
  coeff_t Jhalf = J / 2.;
  coeff_t Jhalf_conj = conj(Jhalf);

  int64_t s1 = op[0];
  int64_t s2 = op[1];

  // Prepare bitmasks
  bit_t flipmask = ((bit_t)1 << s1) | ((bit_t)1 << s2);
  int64_t l = std::min(s1, s2);
  int64_t u = std::max(s1, s2);
  bit_t fermimask = (((bit_t)1 << (u - l - 1)) - 1) << (l + 1);

  auto non_zero_term_ups = [&](bit_t up) {
    return bits::popcnt(up & flipmask) == 1;
  };
  auto non_zero_term_dns = [&](bit_t dn) {
    return bits::popcnt(dn & flipmask) == 1;
  };
  auto term_actionups = [&](bit_t up) -> bit_t { return up ^ flipmask; };
  auto term_actiondns = [&](bit_t dn) -> bit_t { return dn ^ flipmask; };
  auto matrix_element = [&](bit_t up, bit_t dn) -> coeff_t {
    bool fermi_up = bits::popcnt(up & fermimask) & 1;
    bool fermi_dn = bits::popcnt(dn & fermimask) & 1;
    if constexpr (isreal<coeff_t>()) {
      return (fermi_up ^ fermi_dn) ? Jhalf : -Jhalf;
    } else {
      if (bits::gbit(up, s1)) {
        return (fermi_up ^ fermi_dn) ? Jhalf : -Jhalf;
      } else {
        return (fermi_up ^ fermi_dn) ? Jhalf_conj : -Jhalf_conj;
      }
    }
  };

  generic_term_mixed<symmetric, coeff_t>(basis_in, basis_out, non_zero_term_ups,
                                         non_zero_term_dns, term_actionups,
                                         term_actiondns, matrix_element, fill);
}

} // namespace xdiag::basis::tj
