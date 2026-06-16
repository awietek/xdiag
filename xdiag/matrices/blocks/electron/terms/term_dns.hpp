// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>
#include <utility>
#ifdef _OPENMP
#include <omp.h>
#endif

#include <vector>

#include <xdiag/basis/basis_electron.hpp>
#include <xdiag/basis/basis_electron_symmetric.hpp>
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
#else
    auto [begin_up, end_up, idx_up] = utils::thread_range(basis_up, 0, 1);
#endif
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
            auto [dns_out, coeff] = term_action(dns_in);
            int64_t idx_dn_out = basis_dn_out.index(dns_out);
            XDIAG_FILL(base_in + idx_dn_in, base_out + idx_dn_out, -coeff);
          }
          ++idx_dn_in;
        }
      } else {
        int64_t idx_dn_in = 0;
        for (bit_t dns_in : basis_dn_in) {
          if (non_zero_term(dns_in)) {
            auto [dns_out, coeff] = term_action(dns_in);
            int64_t idx_dn_out = basis_dn_out.index(dns_out);
            XDIAG_FILL(base_in + idx_dn_in, base_out + idx_dn_out, coeff);
          }
          ++idx_dn_in;
        }
      }
    }
#ifdef _OPENMP
  }
#endif
}

// Symmetric overload. The dn sector moves while the up sector is a spectator
// representative (basis_in == basis_out: a dn operator conserves Nup and leaves
// the up rep fixed). The dn output must be re-symmetrised within the SAME up
// representative's stabilizer block. The (-1)^Nup Jordan-Wigner factor for an
// odd number of dn operators (fermi_ups) is derived from first principles as
// (-1)^popcount(ups_rep), matching the non-symmetric overload's
// negate = popcount(ups)&1, and is combined with the symmetrisation fermi sign
// (fermi_up XOR fermi_dn) on the non-trivial-stabilizer path.
template <typename coeff_t, bool fermi_ups, typename enumeration_t,
          typename non_zero_term_f, typename term_action_f, typename fill_f>
void term_dns(basis::BasisElectronSymmetric<enumeration_t> const &basis_in,
              basis::BasisElectronSymmetric<enumeration_t> const &basis_out,
              non_zero_term_f non_zero_term, term_action_f term_action,
              fill_f fill) {
  using bit_t = typename enumeration_t::bit_t;
  (void)basis_out; // == basis_in for a dn-sector (Nup-conserving) operator

  auto const &basis_up = basis_in.basis_up();
  auto bloch = basis_in.characters().template as<arma::Col<coeff_t>>();

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
      bool nup_neg = false;
      if constexpr (fermi_ups) {
        nup_neg = bits::popcount(ups) & 1;
      }
      int64_t off = basis_in.ups_offset(idx_up);
      auto dnss = basis_in.dns_for_ups_rep(idx_up);
      auto norms = basis_in.norms_for_ups_rep(idx_up);

      if (basis_in.stab_size(idx_up) == 1) {
        int64_t dn_idx = 0;
        for (bit_t dns_in : dnss) {
          if (non_zero_term(dns_in)) {
            auto [dns_flip, coeff] = term_action(dns_in);
            int64_t idx_dns_flip = basis_in.basis_dn().index(dns_flip);
            XDIAG_FILL(off + dn_idx, off + idx_dns_flip,
                       nup_neg ? -coeff : coeff);
          }
          ++dn_idx;
        }
      } else {
        std::vector<int64_t> syms = basis_in.syms_ups(ups);
        int64_t dn_idx = 0;
        for (bit_t dns_in : dnss) {
          if (non_zero_term(dns_in)) {
            auto [dns_flip, coeff] = term_action(dns_in);
            auto [idx_dns_flip, fermi_dn, sym] =
                basis_in.index_dns_fermi_sym(dns_flip, syms, dnss);
            if (idx_dns_flip >= 0) {
              bool fermi_up = basis_in.fermi_bool_ups(sym, ups);
              coeff_t val = coeff * bloch(sym) * norms[idx_dns_flip] /
                            norms[dn_idx];
              bool neg = (fermi_up ^ fermi_dn) ^ nup_neg;
              XDIAG_FILL(off + dn_idx, off + idx_dns_flip, neg ? -val : val);
            }
          }
          ++dn_idx;
        }
      }
    }
#ifdef _OPENMP
  }
#endif
}

} // namespace xdiag::matrices::electron
