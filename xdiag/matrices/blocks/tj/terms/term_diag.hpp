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
#include <xdiag/matrices/fill_functions.hpp>
#include <xdiag/utils/thread_range.hpp>

namespace xdiag::matrices::tj {

// Compressed rank of physical site s within the non-up complement of `ups`, or
// -1 if s carries an up (so it can hold no dn). Test a dn at s as
// bits::get(dnc, rank). Computed once per ups (hoisted) by the diagonal kernels.
template <typename bit_t>
inline int64_t tj_dn_rank(bit_t ups, int64_t s, int64_t nsites,
                          bit_t sitesmask) {
  if (bits::get(ups, s)) {
    return -1;
  }
  return bits::popcount((~ups) & sitesmask & bits::bitmask<bit_t>(nsites, s));
}

// Diagonal term on the (non-symmetric) tJ basis. A diagonal operator conserves
// both Nup and Ndn, so it only connects a block to itself: when basis_in !=
// basis_out (different number sectors) there is no matrix element.
//
// The diagonal value depends on per-site occupancies, but the dn sector is
// stored compressed, so a physical site s maps to compressed rank
// popcount((~ups) & below(s)). Rather than recompute that (or decompress the dn
// string) for every state, `build(ups)` is called ONCE per ups and returns a
// per-ups evaluator dnc -> coeff that closes over the hoisted ranks; the inner
// loop is then a handful of compressed bit-tests, with no scan on the hot path.
template <typename coeff_t, typename enumeration_t, typename build_f,
          typename fill_f>
void term_diag(basis::BasistJ<enumeration_t> const &basis_in,
               basis::BasistJ<enumeration_t> const &basis_out, build_f build,
               fill_f fill) {
  using bit_t = typename enumeration_t::bit_t;

  if (basis_in != basis_out) {
    return;
  }

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
      auto value = build(ups); // per-ups evaluator: dnc -> coeff (ranks hoisted)
      int64_t idx = basis_in.ups_offset(idx_up);
      for (bit_t dnc : basis_in.basis_dncs(ups)) {
        XDIAG_FILL(idx, idx, value(dnc));
        ++idx;
      }
    }
#ifdef _OPENMP
  }
#endif
}

} // namespace xdiag::matrices::tj
