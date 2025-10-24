// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/common.hpp>

namespace xdiag::basis::electron {

template <typename coeff_t, class basis_t, class apply_f, class fill_f>
void generic_term_sym(basis_t const &basis_in, basis_t const &basis_out,
                      apply_f apply, fill_f fill) {
  using bit_t = typename basis_t::bit_t;

  Representation const &irrep = basis_out.irrep();
  auto characters = irrep.characters().as<arma::Col<coeff_t>>();

  // Loop over all up configurations
#ifdef _OPENMP
#pragma omp parallel
  {
    int num_thread = omp_get_thread_num();
#pragma omp for schedule(runtime)
#endif

    for (int64_t idx_ups = 0; idx_ups < basis_in.n_rep_ups(); ++idx_ups) {
      bit_t ups = basis_in.rep_ups(idx_ups);

      // get the offset of the ups, the corresponding dnss and their norms
      int64_t up_offset_in = basis_in.ups_offset(idx_ups);
      gsl::span<bit_t const> dnss_in = basis_in.dns_for_ups_rep(ups);
      gsl::span<double const> norms_in = basis_in.norms_for_ups_rep(ups);

      int64_t idx_dn = 0;
      for (bit_t dns : dnss_in) {
        // apply the term -> coefficient c and flipped ups/dns
        auto [c, ups_flip, dns_flip] = apply(ups, dns);

        // get the index of the flipped ups, the representative, and the offset
        // of the ups
        int64_t idx_ups_flip = basis_out.index_ups(ups_flip);
        bit_t ups_flip_rep = basis_out.rep_ups(idx_ups_flip);
        int64_t up_offset_out = basis_out.ups_offset(idx_ups_flip);

        // get the symmetries giving the representative, the dns for the rep and
        // their norms
        gsl::span<int64_t const> syms_up_out = basis_out.syms_ups(ups_flip);
        gsl::span<bit_t const> dnss_out =
            basis_out.dns_for_ups_rep(ups_flip_rep);

        // trivial up-stabilizer (likely)
        if (syms_up_out.size() == 1) {
          int64_t sym = syms_up_out.front();

          // get fermi signs from permutation and index of flipped dns
          bool fermi_up = basis_out.fermi_bool_ups(sym, ups_flip);
          auto [idx_dn_flip, fermi_dn] =
              basis_out.index_dns_fermi(dns_flip, sym);

          // out norm is 1 since we have trivial stabilizer
          coeff_t val = c * characters(sym) / norms_in[idx_dn];
          int64_t idx_in = up_offset_in + idx_dn;
          int64_t idx_out = up_offset_out + idx_dn_flip;
          XDIAG_FILL(idx_in, idx_out, (fermi_up ^ fermi_dn) ? -val : val);

        } else { // non-trivial stabilizer
          auto const &syms = syms_up_out;

          gsl::span<double const> norms_out =
              basis_out.norms_for_ups_rep(ups_flip_rep);

          // get fermi signs from permutation and index of flipped dns
          auto [idx_dn_flip, fermi_dn, sym] =
              basis_out.index_dns_fermi_sym(dns_flip, syms, dnss_out);
          bool fermi_up = basis_out.fermi_bool_ups(sym, ups_flip);

          if (idx_dn_flip != invalid_index) {
            coeff_t val =
                c * characters(sym) * norms_out[idx_dn_flip] / norms_in[idx_dn];
            int64_t idx_in = up_offset_in + idx_dn;
            int64_t idx_out = up_offset_out + idx_dn_flip;
            XDIAG_FILL(idx_in, idx_out, (fermi_up ^ fermi_dn) ? -val : val);
          }
        }
        ++idx_dn;
      }
    }
  }
} // namespace xdiag::basis::electron

template <typename coeff_t, class basis_t, class apply_f, class fill_f>
void generic_term_no_sym(basis_t const &basis_in, basis_t const &basis_out,
                         apply_f apply, fill_f fill) {
  using bit_t = typename basis_t::bit_t;
#ifdef _OPENMP
#pragma omp parallel
  {
    int num_thread = omp_get_thread_num();
    auto ups_and_idces = basis_in.states_indices_ups_thread();
#else
  auto ups_and_idces = basis_in.states_indices_ups();
#endif
    int64_t size_dns = basis_out.size_dns();

    for (auto [ups, idx_up] : ups_and_idces) {
      int64_t idx_in = idx_up * size_dns;
      for (auto dns : basis_in.states_dns()) {
        auto [c, ups_flip, dns_flip] = apply(ups, dns);
        int64_t idx_ups_flip = basis_out.index_ups(ups_flip);
        int64_t idx_out_offset = idx_ups_flip * size_dns;
        int64_t idx_out = idx_out_offset + basis_out.index_dns(dns_flip);
        XDIAG_FILL(idx_in, idx_out, c);
      }
      ++idx_in;
#ifdef _OPENMP
    }
#endif
  }
}

template <bool symmetric, typename coeff_t, class basis_t, class apply_f,
          class fill_f>
void generic_term(basis_t const &basis_in, basis_t const &basis_out,
                  apply_f apply, fill_f fill) {
  if constexpr (symmetric) {
    generic_term_sym<coeff_t>(basis_in, basis_out, apply, fill);
  } else {
    generic_term_no_sym<coeff_t>(basis_in, basis_out, apply, fill);
  }
}
} // namespace xdiag::basis::electron
