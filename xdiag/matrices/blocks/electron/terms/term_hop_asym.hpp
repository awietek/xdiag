// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <algorithm>
#include <cstdint>
#include <utility>

#include <xdiag/bits/bitmask.hpp>
#include <xdiag/bits/get_set.hpp>
#include <xdiag/bits/popcount.hpp>
#include <xdiag/bits/zero_one.hpp>
#include <xdiag/matrices/blocks/electron/terms/term_dns.hpp>
#include <xdiag/matrices/blocks/electron/terms/term_ups.hpp>
#include <xdiag/operators/coeff.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/utils/error.hpp>

namespace xdiag::matrices::electron {

// HopupAsym{s1,s2} = -Cdagup{s1} Cup{s2} + Cdagup{s2} Cup{s1}, acting on the up
// sector (matching electron_expansion_rules). The fermi sign maps an expansion
// coefficient -1 to (fermi ? t : -t), as in term_hopup; +1 is its negation.
template <typename coeff_t, class enumeration_t, class fill_f>
void term_hopup_asym(Coeff const &c, Op const &op,
                     basis::BasisElectron<enumeration_t> const &basis_in,
                     basis::BasisElectron<enumeration_t> const &basis_out,
                     fill_f fill) try {
  using bit_t = typename enumeration_t::bit_t;

  int64_t nsites = basis_in.nsites();
  coeff_t t = c.scalar().as<coeff_t>();
  int64_t s1 = op[0];
  int64_t s2 = op[1];
  bit_t flipmask = bits::zero<bit_t>(nsites);
  bits::set(flipmask, s1);
  bits::set(flipmask, s2);
  int64_t l = std::min(s1, s2);
  int64_t u = std::max(s1, s2);
  bit_t fermimask = bits::bitmask<bit_t>(nsites, u - l - 1) << (l + 1);

  auto non_zero_term = [&](bit_t ups) -> bool {
    return bits::popcount(ups & flipmask) & 1;
  };
  auto term_action = [&](bit_t ups) -> std::pair<bit_t, coeff_t> {
    bool fermi = bits::popcount(ups & fermimask) & 1;
    bit_t ups_flip = ups ^ flipmask;
    if (bits::get(ups_flip, s2)) { // Cdagup{s2} Cup{s1}, coefficient +1
      return {ups_flip, fermi ? -t : t};
    } else { // Cdagup{s1} Cup{s2}, coefficient -1
      return {ups_flip, fermi ? t : -t};
    }
  };
  term_ups<coeff_t>(basis_in, basis_out, non_zero_term, term_action, fill);
}
XDIAG_CATCH

// HopdnAsym{s1,s2} = -Cdagdn{s1} Cdn{s2} + Cdagdn{s2} Cdn{s1}, acting on the dn
// sector. Even number of dn operators, so the up parity cancels.
template <typename coeff_t, class enumeration_t, class fill_f>
void term_hopdn_asym(Coeff const &c, Op const &op,
                     basis::BasisElectron<enumeration_t> const &basis_in,
                     basis::BasisElectron<enumeration_t> const &basis_out,
                     fill_f fill) try {
  using bit_t = typename enumeration_t::bit_t;

  int64_t nsites = basis_in.nsites();
  coeff_t t = c.scalar().as<coeff_t>();
  int64_t s1 = op[0];
  int64_t s2 = op[1];
  bit_t flipmask = bits::zero<bit_t>(nsites);
  bits::set(flipmask, s1);
  bits::set(flipmask, s2);
  int64_t l = std::min(s1, s2);
  int64_t u = std::max(s1, s2);
  bit_t fermimask = bits::bitmask<bit_t>(nsites, u - l - 1) << (l + 1);

  auto non_zero_term = [&](bit_t dns) -> bool {
    return bits::popcount(dns & flipmask) & 1;
  };
  auto term_action = [&](bit_t dns) -> std::pair<bit_t, coeff_t> {
    bool fermi = bits::popcount(dns & fermimask) & 1;
    bit_t dns_flip = dns ^ flipmask;
    if (bits::get(dns_flip, s2)) { // Cdagdn{s2} Cdn{s1}, coefficient +1
      return {dns_flip, fermi ? -t : t};
    } else { // Cdagdn{s1} Cdn{s2}, coefficient -1
      return {dns_flip, fermi ? t : -t};
    }
  };
  term_dns<coeff_t, false>(basis_in, basis_out, non_zero_term, term_action,
                           fill);
}
XDIAG_CATCH

} // namespace xdiag::matrices::electron
