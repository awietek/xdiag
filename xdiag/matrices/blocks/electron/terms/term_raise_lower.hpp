// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>
#include <utility>

#include <xdiag/bits/bitmask.hpp>
#include <xdiag/bits/get_set.hpp>
#include <xdiag/bits/nonzero.hpp>
#include <xdiag/bits/popcount.hpp>
#include <xdiag/bits/zero_one.hpp>
#include <xdiag/matrices/blocks/electron/terms/term_dns.hpp>
#include <xdiag/matrices/blocks/electron/terms/term_ups.hpp>
#include <xdiag/operators/coeff.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/utils/error.hpp>

// Single creation / annihilation operators. They change Nup (Cup/Cdagup) or Ndn
// (Cdn/Cdagdn) by one, so basis_in and basis_out are different number sectors.
// The up operators carry the Jordan-Wigner sign of the up occupation below the
// site; the dn operators carry the sign of the dn occupation below the site and,
// because a single dn operator sits behind the entire up string, an additional
// (-1)^Nup (fermi_ups = true).

namespace xdiag::matrices::electron {

template <typename coeff_t, class enumeration_t, class fill_f>
void term_cup(Coeff const &c, Op const &op,
              basis::BasisElectron<enumeration_t> const &basis_in,
              basis::BasisElectron<enumeration_t> const &basis_out,
              fill_f fill) try {
  using bit_t = typename enumeration_t::bit_t;

  int64_t nsites = basis_in.nsites();
  coeff_t cf = c.scalar().as<coeff_t>();
  int64_t s = op[0];
  bit_t mask = bits::zero<bit_t>(nsites);
  bits::set(mask, s);
  bit_t fermi_mask = bits::bitmask<bit_t>(nsites, s);

  auto non_zero_term = [&](bit_t ups) -> bool {
    return bits::nonzero(ups & mask);
  };
  auto term_action = [&](bit_t ups) -> std::pair<bit_t, coeff_t> {
    bool fermi = bits::popcount(ups & fermi_mask) & 1;
    return {ups ^ mask, fermi ? -cf : cf};
  };
  term_ups<coeff_t>(basis_in, basis_out, non_zero_term, term_action, fill);
}
XDIAG_CATCH

template <typename coeff_t, class enumeration_t, class fill_f>
void term_cdagup(Coeff const &c, Op const &op,
                 basis::BasisElectron<enumeration_t> const &basis_in,
                 basis::BasisElectron<enumeration_t> const &basis_out,
                 fill_f fill) try {
  using bit_t = typename enumeration_t::bit_t;

  int64_t nsites = basis_in.nsites();
  coeff_t cf = c.scalar().as<coeff_t>();
  int64_t s = op[0];
  bit_t mask = bits::zero<bit_t>(nsites);
  bits::set(mask, s);
  bit_t fermi_mask = bits::bitmask<bit_t>(nsites, s);

  auto non_zero_term = [&](bit_t ups) -> bool {
    return !bits::nonzero(ups & mask);
  };
  auto term_action = [&](bit_t ups) -> std::pair<bit_t, coeff_t> {
    bool fermi = bits::popcount(ups & fermi_mask) & 1;
    return {ups | mask, fermi ? -cf : cf};
  };
  term_ups<coeff_t>(basis_in, basis_out, non_zero_term, term_action, fill);
}
XDIAG_CATCH

template <typename coeff_t, class enumeration_t, class fill_f>
void term_cdn(Coeff const &c, Op const &op,
              basis::BasisElectron<enumeration_t> const &basis_in,
              basis::BasisElectron<enumeration_t> const &basis_out,
              fill_f fill) try {
  using bit_t = typename enumeration_t::bit_t;

  int64_t nsites = basis_in.nsites();
  coeff_t cf = c.scalar().as<coeff_t>();
  int64_t s = op[0];
  bit_t mask = bits::zero<bit_t>(nsites);
  bits::set(mask, s);
  bit_t fermi_mask = bits::bitmask<bit_t>(nsites, s);

  auto non_zero_term = [&](bit_t dns) -> bool {
    return bits::nonzero(dns & mask);
  };
  auto term_action = [&](bit_t dns) -> std::pair<bit_t, coeff_t> {
    bool fermi = bits::popcount(dns & fermi_mask) & 1;
    return {dns ^ mask, fermi ? -cf : cf};
  };
  term_dns<coeff_t, true>(basis_in, basis_out, non_zero_term, term_action, fill);
}
XDIAG_CATCH

template <typename coeff_t, class enumeration_t, class fill_f>
void term_cdagdn(Coeff const &c, Op const &op,
                 basis::BasisElectron<enumeration_t> const &basis_in,
                 basis::BasisElectron<enumeration_t> const &basis_out,
                 fill_f fill) try {
  using bit_t = typename enumeration_t::bit_t;

  int64_t nsites = basis_in.nsites();
  coeff_t cf = c.scalar().as<coeff_t>();
  int64_t s = op[0];
  bit_t mask = bits::zero<bit_t>(nsites);
  bits::set(mask, s);
  bit_t fermi_mask = bits::bitmask<bit_t>(nsites, s);

  auto non_zero_term = [&](bit_t dns) -> bool {
    return !bits::nonzero(dns & mask);
  };
  auto term_action = [&](bit_t dns) -> std::pair<bit_t, coeff_t> {
    bool fermi = bits::popcount(dns & fermi_mask) & 1;
    return {dns | mask, fermi ? -cf : cf};
  };
  term_dns<coeff_t, true>(basis_in, basis_out, non_zero_term, term_action, fill);
}
XDIAG_CATCH

} // namespace xdiag::matrices::electron
