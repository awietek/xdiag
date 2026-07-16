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
#include <xdiag/kernels/blocks/tj/terms/recompress.hpp>
#include <xdiag/kernels/fill_functions.hpp>
#include <xdiag/operators/coeff.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/thread_range.hpp>

namespace xdiag::kernels::tj {

// Hopup{s1,s2} = -(Cdagup{s1} Cup{s2} + Cdagup{s2} Cup{s1}) on the tJ basis.
//
// Unlike a dn hop, an up hop changes ups, which (a) is only allowed if the
// destination site carries no dn (otherwise the result is doubly occupied and
// not in the basis) and (b) moves the destination out of the dn complement and
// the source into it -- so the compressed dn string must be re-mapped. Because
// both source (had an up) and destination (must be dn-empty) carry no dn, that
// re-mapping is exactly one removed slot + one inserted empty slot, all O(1) and
// free of pdep/pext.
//
// The up Jordan-Wigner sign (popcount(ups & between) & 1) depends on ups only, so
// it is hoisted out of the dn loop; an up hop is two up operators (even), so
// there is no (-1)^Nup cross sign on the dn string.
template <typename coeff_t, class enumeration_t, class fill_f>
void term_hopup(Coeff const &c, Op const &op,
                basis::BasistJ<enumeration_t> const &basis_in,
                basis::BasistJ<enumeration_t> const &basis_out,
                fill_f fill) try {
  using bit_t = typename enumeration_t::bit_t;

  int64_t nsites = basis_in.nsites();
  coeff_t t = c.scalar().as<coeff_t>();
  int64_t s1 = op[0];
  int64_t s2 = op[1];
  bit_t sitesmask = bits::bitmask<bit_t>(nsites, nsites);
  bit_t flipmask = bits::zero<bit_t>(nsites);
  bits::set(flipmask, s1);
  bits::set(flipmask, s2);
  int64_t l = std::min(s1, s2);
  int64_t u = std::max(s1, s2);
  bit_t fermimask = bits::bitmask<bit_t>(nsites, u - l - 1) << (l + 1);

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
      bool up1 = bits::get(ups, s1);
      bool up2 = bits::get(ups, s2);
      if (up1 == up2) {
        continue; // need exactly one up among {s1,s2} to hop
      }

      // --- per-ups setup (hoisted) ---
      int64_t source = up1 ? s1 : s2; // the site that has the up
      int64_t dest = up1 ? s2 : s1;   // the site the up moves to
      bit_t ups_out = ups ^ flipmask;
      bit_t notups_in = (~ups) & sitesmask;
      bit_t notups_out = (~ups_out) & sitesmask;
      // destination rank in the IN complement, source rank in the OUT complement
      int64_t rd_in =
          bits::popcount(notups_in & bits::bitmask<bit_t>(nsites, dest));
      int64_t rs_out =
          bits::popcount(notups_out & bits::bitmask<bit_t>(nsites, source));
      bool fermi = bits::popcount(ups & fermimask) & 1;
      coeff_t coeff = fermi ? t : -t;
      int64_t base_in = basis_in.ups_offset(idx_up);
      int64_t base_out = basis_out.ups_offset(basis_out.index_up(ups_out));

      // --- inner loop: dn is a spectator, but must be dn-empty at dest and
      //     re-mapped between the two compressed spaces ---
      int64_t idx_dnc = 0;
      for (bit_t dnc : basis_in.basis_dncs(ups)) {
        if (!bits::get(dnc, rd_in)) { // destination free of a dn
          bit_t dnc_out =
              insert_zero_slot(remove_slot(dnc, rd_in, nsites), rs_out, nsites);
          int64_t idx_out = base_out + basis_out.index_dncs(ups_out, dnc_out);
          XDIAG_FILL(base_in + idx_dnc, idx_out, coeff);
        }
        ++idx_dnc;
      }
    }
#ifdef _OPENMP
  }
#endif
}
XDIAG_CATCH

} // namespace xdiag::kernels::tj
