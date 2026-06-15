// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>
#include <utility>
#ifdef _OPENMP
#include <omp.h>
#endif

#include <xdiag/basis/basis_electron.hpp>
#include <xdiag/matrices/fill_functions.hpp>
#include <xdiag/utils/thread_range.hpp>

namespace xdiag::matrices::electron {

// Off-diagonal term acting on the up sector only (the dn configuration is left
// untouched). `term_action(ups)` returns {ups_out, coeff}, where coeff already
// carries the fermionic sign of the up-sector action (up operators come first
// in the Jordan-Wigner string, so their sign involves up occupation only).
//
// The dn sector is unchanged, so its enumeration is identical in basis_in and
// basis_out (Ndn is conserved); the dn index is the same offset on both sides.
// The linear index is idx_up * size_dn + idx_dn.
template <typename coeff_t, typename enumeration_t, typename non_zero_term_f,
          typename term_action_f, typename fill_f>
void term_ups(basis::BasisElectron<enumeration_t> const &basis_in,
              basis::BasisElectron<enumeration_t> const &basis_out,
              non_zero_term_f non_zero_term, term_action_f term_action,
              fill_f fill) {
  using bit_t = typename enumeration_t::bit_t;

  auto const &basis_up_in = basis_in.basis_up();
  auto const &basis_up_out = basis_out.basis_up();
  int64_t size_dn_in = basis_in.basis_dn().size();
  int64_t size_dn_out = basis_out.basis_dn().size();

#ifdef _OPENMP
#pragma omp parallel
  {
    int num_thread = omp_get_thread_num();
    auto [begin_up, end_up, idx_up_in] =
        utils::thread_range(basis_up_in, num_thread, omp_get_num_threads());
    for (auto it_up = begin_up; it_up != end_up; ++it_up, ++idx_up_in) {
      bit_t ups_in = *it_up;
      if (non_zero_term(ups_in)) {
        std::pair<bit_t, coeff_t> action = term_action(ups_in);
        bit_t ups_out = action.first;
        coeff_t coeff = action.second;
        int64_t idx_up_out = basis_up_out.index(ups_out);
        int64_t idx_in = idx_up_in * size_dn_in;
        int64_t idx_out = idx_up_out * size_dn_out;
        for (int64_t k = 0; k < size_dn_in; ++k) {
          XDIAG_FILL(idx_in + k, idx_out + k, coeff);
        }
      }
    }
  }
#else
  int64_t idx_up_in = 0;
  for (bit_t ups_in : basis_up_in) {
    if (non_zero_term(ups_in)) {
      std::pair<bit_t, coeff_t> action = term_action(ups_in);
      bit_t ups_out = action.first;
      coeff_t coeff = action.second;
      int64_t idx_up_out = basis_up_out.index(ups_out);
      int64_t idx_in = idx_up_in * size_dn_in;
      int64_t idx_out = idx_up_out * size_dn_out;
      for (int64_t k = 0; k < size_dn_in; ++k) {
        fill(idx_in + k, idx_out + k, coeff);
      }
    }
    ++idx_up_in;
  }
#endif
}

} // namespace xdiag::matrices::electron
