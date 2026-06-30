// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/bits/bitmask.hpp>
#include <xdiag/bits/popcount.hpp>
#include <xdiag/math/complex.hpp>
#include <xdiag/operators/coeff.hpp>
#include <xdiag/operators/op.hpp>

#include <xdiag/kernels/blocks/distributed/electron_distributed/terms/generic_term_dns.hpp>
#include <xdiag/kernels/blocks/distributed/electron_distributed/terms/generic_term_ups.hpp>

namespace xdiag::basis::electron_distributed {

template <typename coeff_t, class basis_t>
void apply_hopping(Coeff const &cpl, Op const &op, basis_t const &basis,
                   const coeff_t *vec_in, coeff_t *vec_out) {
  using bit_t = typename basis_t::bit_t;

  coeff_t t = cpl.scalar().as<coeff_t>();
  std::string type = op.type();
  bool is_up = (type == "Hopup") || (type == "HopupAsym");
  bool is_asym = (type == "HopupAsym") || (type == "HopdnAsym");
  int64_t s1 = op[0];
  int64_t s2 = op[1];

  bit_t flipmask = ((bit_t)1 << s1) | ((bit_t)1 << s2);
  int64_t l = std::min(s1, s2);
  int64_t u = std::max(s1, s2);
  bit_t fermimask = (((bit_t)1 << (u - l - 1)) - 1) << (l + 1);

  // Dispatch on the spin species once; the term_action lambda is selected once
  // below, so the hot loop carries no per-element symmetric/asymmetric branch.
  auto run = [&](auto term_action) {
    if (is_up) {
      auto non_zero_term_dns = [](bit_t dns) -> bool { return true; };
      auto non_zero_term_ups = [&](bit_t ups) -> bool {
        return bits::popcount(ups & flipmask) & 1;
      };
      electron_distributed::generic_term_ups<coeff_t>(
          basis, basis, non_zero_term_ups, non_zero_term_dns, term_action,
          vec_in, vec_out);
    } else {
      auto non_zero_term_ups = [](bit_t ups) -> bool { return true; };
      auto non_zero_term_dns = [&flipmask](bit_t dns) -> bool {
        return bits::popcount(dns & flipmask) & 1;
      };
      electron_distributed::generic_term_dns<coeff_t, false>(
          basis, basis, non_zero_term_ups, non_zero_term_dns, term_action,
          vec_in, vec_out);
    }
  };

  // No conjugation of the coefficient (matching the non-distributed kernels in
  // kernels/blocks/electron). Plain Hop{up,dn} applies the same t in both hop
  // directions; Hop{up,dn}Asym = -Cdag{s1}C{s2} + Cdag{s2}C{s1} flips the sign
  // of the post-flip s2->s1 branch (mirrors electron::term_hop*_asym).
  if (is_asym) {
    run([&](bit_t spins) -> std::pair<bit_t, coeff_t> {
      bool fermi = bits::popcount(spins & fermimask) & 1;
      spins ^= flipmask;
      coeff_t base = fermi ? t : -t;
      return {spins, bits::get(spins, s2) ? -base : base};
    });
  } else {
    run([&](bit_t spins) -> std::pair<bit_t, coeff_t> {
      bool fermi = bits::popcount(spins & fermimask) & 1;
      spins ^= flipmask;
      return {spins, fermi ? t : -t};
    });
  }
}

} // namespace xdiag::basis::electron_distributed
