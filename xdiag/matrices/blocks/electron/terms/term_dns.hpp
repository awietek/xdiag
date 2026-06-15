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
#include <xdiag/bits/popcount.hpp>
#include <xdiag/matrices/fill_functions.hpp>
#include <xdiag/utils/thread_range.hpp>

namespace xdiag::matrices::electron {

// Off-diagonal term acting on the dn sector only (the up configuration is left
// untouched). `term_action(dns)` returns {dns_out, coeff}, where coeff carries
// the fermionic sign of the dn-sector action (dn occupation between/below the
// involved sites).
//
// In the Jordan-Wigner ordering "all ups then all dns", a dn operator string is
// preceded by the entire up string. For an action with an ODD number of dn
// fermion operators (a single Cdn / Cdagdn) this contributes an extra sign
// (-1)^(number of up electrons); set `fermi_ups = true` for those. For an EVEN
// number (a dn hopping Cdagdn Cdn) the up parity cancels; use `fermi_ups =
// false`.
//
// The outer loop runs over the up sector (a spectator, Nup conserved, identical
// enumeration in/out). The up-parity sign is decided once per up state and the
// inner dn loop is selected accordingly, so no sign work happens in the hot
// loop. For fermi_ups == false the `negate` branch folds away entirely and only
// the +coeff loop remains. The linear index is idx_up * size_dn + idx_dn.
template <typename coeff_t, bool fermi_ups, typename enumeration_t,
          typename non_zero_term_f, typename term_action_f, typename fill_f>
void term_dns(basis::BasisElectron<enumeration_t> const &basis_in,
              basis::BasisElectron<enumeration_t> const &basis_out,
              non_zero_term_f non_zero_term, term_action_f term_action,
              fill_f fill) {
  using bit_t = typename enumeration_t::bit_t;

  auto const &basis_up = basis_in.basis_up();
  auto const &basis_dn_in = basis_in.basis_dn();
  auto const &basis_dn_out = basis_out.basis_dn();
  int64_t size_dn_in = basis_dn_in.size();
  int64_t size_dn_out = basis_dn_out.size();

#ifdef _OPENMP
#pragma omp parallel
  {
    int num_thread = omp_get_thread_num();
    auto [begin_up, end_up, idx_up] =
        utils::thread_range(basis_up, num_thread, omp_get_num_threads());
    for (auto it_up = begin_up; it_up != end_up; ++it_up, ++idx_up) {
      bit_t ups = *it_up;
      int64_t base_in = idx_up * size_dn_in;
      int64_t base_out = idx_up * size_dn_out;
      bool negate = false;
      if constexpr (fermi_ups) {
        negate = bits::popcount(ups) & 1;
      }
      if (negate) {
        int64_t idx_dn_in = 0;
        for (bit_t dns_in : basis_dn_in) {
          if (non_zero_term(dns_in)) {
            std::pair<bit_t, coeff_t> action = term_action(dns_in);
            int64_t idx_dn_out = basis_dn_out.index(action.first);
            XDIAG_FILL(base_in + idx_dn_in, base_out + idx_dn_out,
                       -action.second);
          }
          ++idx_dn_in;
        }
      } else {
        int64_t idx_dn_in = 0;
        for (bit_t dns_in : basis_dn_in) {
          if (non_zero_term(dns_in)) {
            std::pair<bit_t, coeff_t> action = term_action(dns_in);
            int64_t idx_dn_out = basis_dn_out.index(action.first);
            XDIAG_FILL(base_in + idx_dn_in, base_out + idx_dn_out,
                       action.second);
          }
          ++idx_dn_in;
        }
      }
    }
  }
#else
  int64_t idx_up = 0;
  for (bit_t ups : basis_up) {
    int64_t base_in = idx_up * size_dn_in;
    int64_t base_out = idx_up * size_dn_out;
    bool negate = false;
    if constexpr (fermi_ups) {
      negate = bits::popcount(ups) & 1;
    }
    if (negate) {
      int64_t idx_dn_in = 0;
      for (bit_t dns_in : basis_dn_in) {
        if (non_zero_term(dns_in)) {
          std::pair<bit_t, coeff_t> action = term_action(dns_in);
          int64_t idx_dn_out = basis_dn_out.index(action.first);
          fill(base_in + idx_dn_in, base_out + idx_dn_out, -action.second);
        }
        ++idx_dn_in;
      }
    } else {
      int64_t idx_dn_in = 0;
      for (bit_t dns_in : basis_dn_in) {
        if (non_zero_term(dns_in)) {
          std::pair<bit_t, coeff_t> action = term_action(dns_in);
          int64_t idx_dn_out = basis_dn_out.index(action.first);
          fill(base_in + idx_dn_in, base_out + idx_dn_out, action.second);
        }
        ++idx_dn_in;
      }
    }
    ++idx_up;
  }
#endif
}

} // namespace xdiag::matrices::electron
