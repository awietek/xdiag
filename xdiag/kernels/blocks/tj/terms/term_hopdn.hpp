// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <algorithm>
#include <cstdint>
#ifdef _OPENMP
#include <omp.h>
#endif

#include <xdiag/basis/basis_tj.hpp>
#include <xdiag/bits/bitmask.hpp>
#include <xdiag/bits/get_set.hpp>
#include <xdiag/bits/popcount.hpp>
#include <xdiag/bits/zero_one.hpp>
#include <xdiag/kernels/fill_functions.hpp>
#include <xdiag/operators/coeff.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/thread_range.hpp>

namespace xdiag::kernels::tj {

// Hopdn{s1,s2} = -(Cdagdn{s1} Cdn{s2} + Cdagdn{s2} Cdn{s1}) on the tJ basis.
//
// The dn sector is stored compressed into the nsites-nup non-up sites, and a dn
// hop conserves ups -- so for each ups this is just a spinless-fermion hop on
// the compressed dn string. The only ups-dependent work is mapping the two
// physical sites to their compressed ranks (two popcounts, hoisted out of the
// dn loop); after that the inner loop is pure compressed-space bit arithmetic,
// with no pdep/pext.
//
// A dn can neither sit on nor hop through a site occupied by an up, so if either
// hop site carries an up the bond contributes nothing for that ups.
template <typename coeff_t, class enumeration_t, class fill_f>
void term_hopdn(Coeff const &c, Op const &op,
                basis::BasistJ<enumeration_t> const &basis_in,
                basis::BasistJ<enumeration_t> const &basis_out,
                fill_f fill) try {
  using bit_t = typename enumeration_t::bit_t;

  int64_t nsites = basis_in.nsites();
  coeff_t t = c.scalar().as<coeff_t>();
  int64_t s1 = op[0];
  int64_t s2 = op[1];
  bit_t sitesmask = bits::bitmask<bit_t>(nsites, nsites);

  auto const &basis_up = basis_in.basis_up();

#ifdef _OPENMP
#pragma omp parallel
  {
    int num_thread = omp_get_thread_num();
    auto [begin_up, end_up, idx_up] =
        utils::thread_range(basis_up, num_thread, omp_get_num_threads());
#else
    auto [begin_up, end_up, idx_up] = utils::thread_range(basis_up, 0, 1);
#endif
    for (auto it_up = begin_up; it_up != end_up; ++it_up, ++idx_up) {
      bit_t ups = *it_up;
      if (bits::get(ups, s1) || bits::get(ups, s2)) {
        continue; // a hop site is up-occupied -> no dn contribution
      }

      // --- per-ups setup (hoisted): physical sites -> compressed ranks ---
      bit_t notups = (~ups) & sitesmask;
      int64_t r1 = bits::popcount(notups & bits::bitmask<bit_t>(nsites, s1));
      int64_t r2 = bits::popcount(notups & bits::bitmask<bit_t>(nsites, s2));
      int64_t rl = std::min(r1, r2);
      int64_t rh = std::max(r1, r2);
      bit_t flipmask = bits::zero<bit_t>(nsites);
      bits::set(flipmask, r1);
      bits::set(flipmask, r2);
      bit_t fermimask = bits::bitmask<bit_t>(nsites, rh - rl - 1) << (rl + 1);

      // --- inner loop: spinless-fermion hop on the compressed dn string ---
      int64_t base = basis_in.ups_offset(idx_up); // ups unchanged by a dn hop
      int64_t idx_dnc_in = 0;
      for (bit_t dnc : basis_in.basis_dncs(ups)) {
        if (bits::popcount(dnc & flipmask) & 1) { // exactly one end occupied
          bit_t dnc_out = dnc ^ flipmask;
          bool fermi = bits::popcount(dnc & fermimask) & 1;
          int64_t idx_dnc_out = basis_out.index_dncs(ups, dnc_out);
          XDIAG_FILL(base + idx_dnc_in, base + idx_dnc_out, fermi ? t : -t);
        }
        ++idx_dnc_in;
      }
    }
#ifdef _OPENMP
  }
#endif
}
XDIAG_CATCH

} // namespace xdiag::kernels::tj
