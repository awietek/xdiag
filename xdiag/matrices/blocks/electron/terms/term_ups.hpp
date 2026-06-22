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
#else
  auto [begin_up, end_up, idx_up_in] = utils::thread_range(basis_up_in, 0, 1);
#endif
    for (auto it_up = begin_up; it_up != end_up; ++it_up, ++idx_up_in) {
      bit_t ups_in = *it_up;
      if (non_zero_term(ups_in)) {
        auto [ups_out, coeff] = term_action(ups_in);
        int64_t idx_up_out = basis_up_out.index(ups_out);
        int64_t idx_in = idx_up_in * size_dn_in;
        int64_t idx_out = idx_up_out * size_dn_out;
        for (int64_t k = 0; k < size_dn_in; ++k) {
          XDIAG_FILL(idx_in + k, idx_out + k, coeff);
        }
      }
    }
#ifdef _OPENMP
  }
#endif
}

// Symmetric overload. The up sector moves (ups_in -> ups_out) while the dn
// configuration is a spectator, but under symmetry the dn must be
// re-symmetrised in the OUTPUT up representative's block. Outer loop over input
// up representatives; the output up representative and its up-stabilizer are
// found once per input up (hoisted out of the inner dn loop). The matrix
// element carries the character of the mapping symmetry, the norm ratio, and
// the fermi sign fermi_up XOR fermi_dn. `term_action`'s coeff already carries
// the up-sector operator sign.
template <typename coeff_t, typename basis_t, typename non_zero_term_f,
          typename term_action_f, typename fill_f,
          typename = decltype(std::declval<basis_t>().dns_for_ups_rep(0))>
void term_ups(basis_t const &basis_in, basis_t const &basis_out,
              non_zero_term_f non_zero_term, term_action_f term_action,
              fill_f fill) {
  using bit_t = typename basis_t::bit_t;

  auto const &basis_up_in = basis_in.basis_up();

  // conjugation necessary for definition of projected states
  arma::Col<coeff_t> characters =
      arma::conj(basis_out.characters().template as<arma::Col<coeff_t>>());

#ifdef _OPENMP
#pragma omp parallel
  {
    int num_thread = omp_get_thread_num();
    auto [begin_up, end_up, idx_up_in] =
        utils::thread_range(basis_up_in, num_thread, omp_get_num_threads());
#else
  auto [begin_up, end_up, idx_up_in] = utils::thread_range(basis_up_in, 0, 1);
#endif
    for (auto it_up = begin_up; it_up != end_up; ++it_up, ++idx_up_in) {
      bit_t ups_in = *it_up;
      if (!non_zero_term(ups_in)) {
        continue;
      }
      auto [ups_out, coeff] = term_action(ups_in);

      auto [raw, s0, nrm] = basis_out.basis_up().representative_data(ups_out);
      (void)nrm;
      int64_t idx_up_out = raw - 1;
      int64_t off_in = basis_in.ups_offset(idx_up_in);
      int64_t off_out = basis_out.ups_offset(idx_up_out);
      auto dnss_in = basis_in.dns_for_ups_rep(idx_up_in);
      auto norms_in = basis_in.norms_for_ups_rep(idx_up_in);

      if (basis_out.stab_size(idx_up_out) == 1) {
        // trivial out-stabilizer: single mapping symmetry s0, norm_out == 1
        coeff_t prefac = coeff * characters(s0);
        bool fermi_up = basis_out.fermi_bool_ups(s0, ups_out);
        auto dnss_out = basis_out.dns_for_ups_rep(idx_up_out);
        int64_t idx_dn = 0;
        for (bit_t dns : dnss_in) {
          auto [idx_dns_out, fermi_dn] =
              basis_out.index_dns_fermi(dns, s0, idx_up_out, dnss_out);
          if (idx_dns_out >=
              0) { // -1: tJ double occupancy (never for electron)
            coeff_t val = prefac / norms_in[idx_dn];
            XDIAG_FILL(off_in + idx_dn, off_out + idx_dns_out,
                       (fermi_up ^ fermi_dn) ? -val : val);
          }
          ++idx_dn;
        }
      } else {
        // non-trivial out-stabilizer
        std::vector<int64_t> syms = basis_out.syms_ups(ups_out);
        auto dnss_out = basis_out.dns_for_ups_rep(idx_up_out);
        auto norms_out = basis_out.norms_for_ups_rep(idx_up_out);
        int64_t idx_dn = 0;
        for (bit_t dns : dnss_in) {
          auto [idx_dns_out, fermi_dn, sym] =
              basis_out.index_dns_fermi_sym(dns, syms, dnss_out);
          if (idx_dns_out >= 0) {
            bool fermi_up = basis_out.fermi_bool_ups(sym, ups_out);
            coeff_t val = coeff * characters(sym) * norms_out[idx_dns_out] /
                          norms_in[idx_dn];
            XDIAG_FILL(off_in + idx_dn, off_out + idx_dns_out,
                       (fermi_up ^ fermi_dn) ? -val : val);
          }
          ++idx_dn;
        }
      }
    }
#ifdef _OPENMP
  }
#endif
}

} // namespace xdiag::matrices::electron
