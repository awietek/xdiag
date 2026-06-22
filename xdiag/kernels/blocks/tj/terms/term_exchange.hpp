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

// Exchange{i,j} = 1/2 (S+{i} S-{j} + S-{i} S+{j}), the off-diagonal spin-flip
// part of the (tJ/electron) Heisenberg coupling SdotS = SzSz + Exchange (the 1/2
// is part of the operator definition, matching the Electron block). It connects
// the antiparallel singly-occupied
// configurations (up at one site, dn at the other) and swaps the two spins; the
// charge configuration is unchanged.
//
// For a given ups exactly one of {i,j} carries an up (else the two sites are not
// antiparallel singly occupied and the term is zero). Call that the up site `a`
// and the other the non-up site `b`; the term acts only if b carries a dn. The
// up change ups -> ups ^ {a,b} is the same as Hopup, so the dn re-mapping reuses
// the Hopup slot surgery -- except the moved slot is filled: the dn leaves b and
// appears at a. The sign is the up hop sign (a->b) times the dn hop sign (b->a)
// times an overall (-1) from the operator ordering of the spin-flip string;
// since the charge config is preserved there is no (-1)^Nup cross sign.
template <typename coeff_t, class enumeration_t, class fill_f>
void term_exchange(Coeff const &c, Op const &op,
                   basis::BasistJ<enumeration_t> const &basis_in,
                   basis::BasistJ<enumeration_t> const &basis_out,
                   fill_f fill) try {
  using bit_t = typename enumeration_t::bit_t;

  int64_t nsites = basis_in.nsites();
  coeff_t J = c.scalar().as<coeff_t>();
  int64_t s1 = op[0];
  int64_t s2 = op[1];
  coeff_t Jhalf = J * coeff_t(0.5); // the 1/2 in 1/2(S+S- + S-S+)
  bit_t sitesmask = bits::bitmask<bit_t>(nsites, nsites);
  bit_t flipmask = bits::zero<bit_t>(nsites);
  bits::set(flipmask, s1);
  bits::set(flipmask, s2);
  int64_t l = std::min(s1, s2);
  int64_t h = std::max(s1, s2);
  bit_t between = bits::bitmask<bit_t>(nsites, h - l - 1) << (l + 1);

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
        continue; // need exactly one up among {s1,s2}
      }

      // --- per-ups setup (hoisted) ---
      int64_t a = up1 ? s1 : s2;    // the up site (loses up, gains dn)
      int64_t b = up1 ? s2 : s1;    // the non-up site (must hold the dn)
      bit_t ups_out = ups ^ flipmask;
      bit_t notups_in = (~ups) & sitesmask;
      bit_t notups_out = (~ups_out) & sitesmask;
      int64_t rb_in =
          bits::popcount(notups_in & bits::bitmask<bit_t>(nsites, b));
      int64_t ra_out =
          bits::popcount(notups_out & bits::bitmask<bit_t>(nsites, a));
      // the dns strictly between a and b are the non-up sites between them; in
      // compressed-in space they occupy `cnt` ranks adjacent to rb_in, on the
      // side toward a -> the dn hop (b->a) Jordan-Wigner mask.
      int64_t cnt = bits::popcount(notups_in & between);
      bit_t dnfmask = (b < a) ? (bits::bitmask<bit_t>(nsites, cnt) << (rb_in + 1))
                              : (bits::bitmask<bit_t>(nsites, cnt)
                                 << (rb_in - cnt));
      bool up_fermi = bits::popcount(ups & between) & 1;
      int64_t base_in = basis_in.ups_offset(idx_up);
      int64_t base_out = basis_out.ups_offset(basis_out.index_up(ups_out));

      int64_t idx_dnc = 0;
      for (bit_t dnc : basis_in.basis_dncs(ups)) {
        if (bits::get(dnc, rb_in)) { // b carries a dn to flip
          bool dn_fermi = bits::popcount(dnc & dnfmask) & 1;
          // remove the dn slot at b, insert a filled dn slot at a
          bit_t dnc_out =
              insert_zero_slot(remove_slot(dnc, rb_in, nsites), ra_out, nsites);
          bits::set(dnc_out, ra_out);
          int64_t idx_out = base_out + basis_out.index_dncs(ups_out, dnc_out);
          bool neg = !(up_fermi ^ dn_fermi); // overall (-1) from spin-flip order
          XDIAG_FILL(base_in + idx_dnc, idx_out, neg ? -Jhalf : Jhalf);
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
