// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>
#ifdef _OPENMP
#include <omp.h>
#endif

#include <xdiag/basis/basis_electron.hpp>
#include <xdiag/matrices/fill_functions.hpp>
#include <xdiag/utils/thread_range.hpp>

namespace xdiag::matrices::electron {

// Diagonal term on the (non-symmetric) electron product basis. `apply(ups, dns)`
// returns the diagonal coefficient for the product state (ups, dns). The linear
// index is idx_up * size_dn + idx_dn (up sector outer, dn sector inner),
// matching BasisElectron::index.
//
// A diagonal operator conserves both Nup and Ndn, so it only connects a block
// to itself: when basis_in != basis_out (different number sectors) there is no
// matrix element and we return immediately.
template <typename coeff_t, typename enumeration_t, typename apply_f,
          typename fill_f>
void term_diag(basis::BasisElectron<enumeration_t> const &basis_in,
               basis::BasisElectron<enumeration_t> const &basis_out,
               apply_f apply, fill_f fill) {
  using bit_t = typename enumeration_t::bit_t;

  if (basis_in != basis_out) {
    return;
  }

  auto const &basis_up = basis_in.basis_up();
  auto const &basis_dn = basis_in.basis_dn();
  int64_t size_dn = basis_dn.size();

#ifdef _OPENMP
#pragma omp parallel
  {
    int num_thread = omp_get_thread_num();
    auto [begin_up, end_up, idx_up] =
        utils::thread_range(basis_up, num_thread, omp_get_num_threads());
    for (auto it_up = begin_up; it_up != end_up; ++it_up, ++idx_up) {
      bit_t ups = *it_up;
      int64_t idx = idx_up * size_dn;
      for (bit_t dns : basis_dn) {
        coeff_t c = apply(ups, dns);
        XDIAG_FILL(idx, idx, c);
        ++idx;
      }
    }
  }
#else
  int64_t idx_up = 0;
  for (bit_t ups : basis_up) {
    int64_t idx = idx_up * size_dn;
    for (bit_t dns : basis_dn) {
      coeff_t c = apply(ups, dns);
      fill(idx, idx, c);
      ++idx;
    }
    ++idx_up;
  }
#endif
}

} // namespace xdiag::matrices::electron
