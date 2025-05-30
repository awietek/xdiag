// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/bits/bitops.hpp>
#include <xdiag/common.hpp>
#include <xdiag/operators/op.hpp>

#include <xdiag/basis/tj_distributed/apply/generic_term_dns.hpp>
#include <xdiag/basis/tj_distributed/apply/generic_term_ups.hpp>

namespace xdiag::basis::tj_distributed {

template <typename coeff_t, class basis_t>
void apply_hopping(Coupling const &cpl, Op const &op, basis_t const &basis,
                   const coeff_t *vec_in, coeff_t *vec_out) {
  using bit_t = typename basis_t::bit_t;

  coeff_t t = cpl.scalar().as<coeff_t>();
  std::string type = op.type();
  int64_t s1 = op[0];
  int64_t s2 = op[1];

  bit_t flipmask = ((bit_t)1 << s1) | ((bit_t)1 << s2);
  int64_t l = std::min(s1, s2);
  int64_t u = std::max(s1, s2);
  bit_t fermimask = (((bit_t)1 << (u - l - 1)) - 1) << (l + 1);

  auto term_action = [&](bit_t spins) -> std::pair<bit_t, coeff_t> {
    bool fermi = bits::popcnt(spins & fermimask) & 1;
    spins ^= flipmask;
    if constexpr (isreal<coeff_t>()) {
      return {spins, fermi ? t : -t};
    } else {
      coeff_t tt = (bits::gbit(spins, s1)) ? t : conj(t);
      return {spins, fermi ? tt : -tt};
    }
  };

  if (type == "Hopup") {

    // Define annihilation conditions
    auto non_zero_term_dns = [&](bit_t const &dns) -> bool {
      return (dns & flipmask) == 0;
    };
    auto non_zero_term_ups = [&](bit_t const &ups) -> bool {
      return bits::popcnt(ups & flipmask) & 1;
    };

    // Call generic term function
    tj_distributed::generic_term_ups<coeff_t>(basis, basis, non_zero_term_ups,
                                              non_zero_term_dns, term_action,
                                              vec_in, vec_out);
  } else if (type == "Hopdn") {

    // Define annihilation conditions
    auto non_zero_term_ups = [&](bit_t const &ups) -> bool {
      return (ups & flipmask) == 0;
    };
    auto non_zero_term_dns = [&](bit_t const &dns) -> bool {
      return bits::popcnt(dns & flipmask) & 1;
    };

    // Call generic term function
    tj_distributed::generic_term_dns<coeff_t, false>(
        basis, basis, non_zero_term_ups, non_zero_term_dns, term_action, vec_in,
        vec_out);
  }
}

} // namespace xdiag::basis::tj_distributed
