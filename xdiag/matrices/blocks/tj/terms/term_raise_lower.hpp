// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>
#ifdef _OPENMP
#include <omp.h>
#endif

#include <xdiag/basis/basis_tj.hpp>
#include <xdiag/bits/bitmask.hpp>
#include <xdiag/bits/get_set.hpp>
#include <xdiag/bits/popcount.hpp>
#include <xdiag/bits/zero_one.hpp>
#include <xdiag/matrices/blocks/tj/terms/recompress.hpp>
#include <xdiag/matrices/fill_functions.hpp>
#include <xdiag/operators/coeff.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/thread_range.hpp>

// Single creation / annihilation operators on the tJ basis. They change Nup
// (Cup/Cdagup) or Ndn (Cdn/Cdagdn) by one, so basis_in and basis_out are
// different number sectors.
//
//   dn operators (Cdn/Cdagdn): ups is unchanged, so the dn complement is
//   unchanged and the compressed dn string just clears/sets the bit at the
//   site's rank -- no re-mapping. They carry the dn Jordan-Wigner sign below the
//   site (in compressed space) and, because a single dn operator sits behind the
//   whole up string, an additional (-1)^Nup.
//
//   up operators (Cup/Cdagup): ups changes, so the site joins (Cup) or leaves
//   (Cdagup) the dn complement and the spectator dn string is re-mapped with one
//   O(1) slot insertion/removal (the slot is always empty: the site carries an
//   up before/after, so no dn). They carry only the up Jordan-Wigner sign below
//   the site. Cdagup additionally requires the site to be dn-empty (else the
//   result is doubly occupied and not in the basis).

namespace xdiag::matrices::tj {

// Cdn{s}: annihilate a dn at site s (Ndn -> Ndn-1).
template <typename coeff_t, class enumeration_t, class fill_f>
void term_cdn(Coeff const &c, Op const &op,
              basis::BasistJ<enumeration_t> const &basis_in,
              basis::BasistJ<enumeration_t> const &basis_out, fill_f fill) try {
  using bit_t = typename enumeration_t::bit_t;
  int64_t nsites = basis_in.nsites();
  coeff_t cf = c.scalar().as<coeff_t>();
  int64_t s = op[0];
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
      if (bits::get(ups, s)) {
        continue; // s carries an up -> no dn there to annihilate
      }
      bit_t notups = (~ups) & sitesmask;
      int64_t r = bits::popcount(notups & bits::bitmask<bit_t>(nsites, s));
      bit_t mask = bits::zero<bit_t>(nsites);
      bits::set(mask, r);
      bit_t fermimask = bits::bitmask<bit_t>(nsites, r); // dn below r
      bool nup_neg = bits::popcount(ups) & 1;
      int64_t base_in = basis_in.ups_offset(idx_up);
      int64_t base_out = basis_out.ups_offset(idx_up); // ups (and Nup) unchanged
      int64_t idx_dnc = 0;
      for (bit_t dnc : basis_in.basis_dncs(ups)) {
        if (bits::nonzero(dnc & mask)) { // dn present at s
          bool sign = (bool)(bits::popcount(dnc & fermimask) & 1) ^ nup_neg;
          int64_t idx_out = base_out + basis_out.index_dncs(ups, dnc ^ mask);
          XDIAG_FILL(base_in + idx_dnc, idx_out, sign ? -cf : cf);
        }
        ++idx_dnc;
      }
    }
#ifdef _OPENMP
  }
#endif
}
XDIAG_CATCH

// Cdagdn{s}: create a dn at site s (Ndn -> Ndn+1).
template <typename coeff_t, class enumeration_t, class fill_f>
void term_cdagdn(Coeff const &c, Op const &op,
                 basis::BasistJ<enumeration_t> const &basis_in,
                 basis::BasistJ<enumeration_t> const &basis_out,
                 fill_f fill) try {
  using bit_t = typename enumeration_t::bit_t;
  int64_t nsites = basis_in.nsites();
  coeff_t cf = c.scalar().as<coeff_t>();
  int64_t s = op[0];
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
      if (bits::get(ups, s)) {
        continue; // s carries an up -> creating a dn would doubly occupy
      }
      bit_t notups = (~ups) & sitesmask;
      int64_t r = bits::popcount(notups & bits::bitmask<bit_t>(nsites, s));
      bit_t mask = bits::zero<bit_t>(nsites);
      bits::set(mask, r);
      bit_t fermimask = bits::bitmask<bit_t>(nsites, r); // dn below r
      bool nup_neg = bits::popcount(ups) & 1;
      int64_t base_in = basis_in.ups_offset(idx_up);
      int64_t base_out = basis_out.ups_offset(idx_up);
      int64_t idx_dnc = 0;
      for (bit_t dnc : basis_in.basis_dncs(ups)) {
        if (bits::iszero(dnc & mask)) { // s free of a dn
          bool sign = (bool)(bits::popcount(dnc & fermimask) & 1) ^ nup_neg;
          int64_t idx_out = base_out + basis_out.index_dncs(ups, dnc | mask);
          XDIAG_FILL(base_in + idx_dnc, idx_out, sign ? -cf : cf);
        }
        ++idx_dnc;
      }
    }
#ifdef _OPENMP
  }
#endif
}
XDIAG_CATCH

// Cup{s}: annihilate an up at site s (Nup -> Nup-1). The site joins the dn
// complement, so the spectator dn string gains an empty slot at its rank.
template <typename coeff_t, class enumeration_t, class fill_f>
void term_cup(Coeff const &c, Op const &op,
              basis::BasistJ<enumeration_t> const &basis_in,
              basis::BasistJ<enumeration_t> const &basis_out, fill_f fill) try {
  using bit_t = typename enumeration_t::bit_t;
  int64_t nsites = basis_in.nsites();
  coeff_t cf = c.scalar().as<coeff_t>();
  int64_t s = op[0];
  bit_t sitesmask = bits::bitmask<bit_t>(nsites, nsites);
  bit_t mask = bits::zero<bit_t>(nsites);
  bits::set(mask, s);
  bit_t fermimask = bits::bitmask<bit_t>(nsites, s); // up below s
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
      if (!bits::get(ups, s)) {
        continue; // no up at s to annihilate
      }
      bit_t ups_out = ups ^ mask;
      bit_t notups_out = (~ups_out) & sitesmask;
      int64_t r_out =
          bits::popcount(notups_out & bits::bitmask<bit_t>(nsites, s));
      bool fermi = bits::popcount(ups & fermimask) & 1;
      coeff_t coeff = fermi ? -cf : cf;
      int64_t base_in = basis_in.ups_offset(idx_up);
      int64_t base_out = basis_out.ups_offset(basis_out.index_up(ups_out));
      int64_t idx_dnc = 0;
      for (bit_t dnc : basis_in.basis_dncs(ups)) {
        bit_t dnc_out = insert_zero_slot(dnc, r_out, nsites);
        int64_t idx_out = base_out + basis_out.index_dncs(ups_out, dnc_out);
        XDIAG_FILL(base_in + idx_dnc, idx_out, coeff);
        ++idx_dnc;
      }
    }
#ifdef _OPENMP
  }
#endif
}
XDIAG_CATCH

// Cdagup{s}: create an up at site s (Nup -> Nup+1). The site leaves the dn
// complement, so the spectator dn string drops its (empty) slot at that rank;
// the site must be dn-empty, else the result is doubly occupied.
template <typename coeff_t, class enumeration_t, class fill_f>
void term_cdagup(Coeff const &c, Op const &op,
                 basis::BasistJ<enumeration_t> const &basis_in,
                 basis::BasistJ<enumeration_t> const &basis_out,
                 fill_f fill) try {
  using bit_t = typename enumeration_t::bit_t;
  int64_t nsites = basis_in.nsites();
  coeff_t cf = c.scalar().as<coeff_t>();
  int64_t s = op[0];
  bit_t sitesmask = bits::bitmask<bit_t>(nsites, nsites);
  bit_t mask = bits::zero<bit_t>(nsites);
  bits::set(mask, s);
  bit_t fermimask = bits::bitmask<bit_t>(nsites, s); // up below s
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
      if (bits::get(ups, s)) {
        continue; // s already carries an up
      }
      bit_t notups = (~ups) & sitesmask;
      int64_t r_in = bits::popcount(notups & bits::bitmask<bit_t>(nsites, s));
      bit_t ups_out = ups | mask;
      bool fermi = bits::popcount(ups & fermimask) & 1;
      coeff_t coeff = fermi ? -cf : cf;
      int64_t base_in = basis_in.ups_offset(idx_up);
      int64_t base_out = basis_out.ups_offset(basis_out.index_up(ups_out));
      int64_t idx_dnc = 0;
      for (bit_t dnc : basis_in.basis_dncs(ups)) {
        if (!bits::get(dnc, r_in)) { // s free of a dn -> result not doubly occ.
          bit_t dnc_out = remove_slot(dnc, r_in, nsites);
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

} // namespace xdiag::matrices::tj
