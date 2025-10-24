// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <functional>
#include <xdiag/parallel/omp/omp_utils.hpp>

namespace xdiag::basis::electron {

template <typename coeff_t, bool fermi_ups, class basis_t,
          class non_zero_term_f, class term_action_f, class fill_f>
void generic_term_dns_sym(basis_t const &basis_in, basis_t const &basis_out,
                          non_zero_term_f non_zero_term,
                          term_action_f term_action, fill_f fill) {
  using bit_t = typename basis_t::bit_t;

  Representation const &irrep = basis_out.irrep();
  auto bloch_factors = irrep.characters().as<arma::Col<coeff_t>>();

#ifdef _OPENMP
#pragma omp parallel
  {
    int num_thread = omp_get_thread_num();
#pragma omp for schedule(runtime)
#endif
    for (int64_t idx_up = 0; idx_up < basis_in.n_rep_ups(); ++idx_up) {
      bit_t ups_in = basis_in.rep_ups(idx_up);
      bit_t ups_out = basis_out.rep_ups(idx_up);
      assert(ups_out == ups_in);

      int64_t ups_offset_in = basis_out.ups_offset(idx_up);
      auto dnss_in = basis_in.dns_for_ups_rep(ups_in);
      auto norms_in = basis_in.norms_for_ups_rep(ups_in);

      int64_t ups_offset_out = basis_out.ups_offset(idx_up);
      auto syms_out = basis_out.syms_ups(ups_out);
      auto dnss_out = basis_out.dns_for_ups_rep(ups_out);
      auto norms_out = basis_out.norms_for_ups_rep(ups_out);

      // trivial up-stabilizer (likely)
      if (syms_out.size() == 1) {
        int64_t dns_in_idx = 0;
        for (bit_t dns_in : dnss_in) {

          if (non_zero_term(dns_in)) {
            int64_t idx_in = ups_offset_in + dns_in_idx;
            auto [dns_flip, coeff] = term_action(dns_in);

            if constexpr (fermi_ups) { // not ideal to do this here
              if (bits::popcnt(ups_in) & 1) {
                coeff = -coeff;
              }
            }

            int64_t idx_dns_flip = basis_out.index_dns(dns_flip);
            int64_t idx_out = ups_offset_out + idx_dns_flip;
            XDIAG_FILL(idx_in, idx_out, coeff);
          }
          ++dns_in_idx;
        }
      }

      else { // non-trivial up-stabilizer (unlikely)

        int64_t dns_in_idx = 0;
        for (bit_t dns_in : dnss_in) {

          if (non_zero_term(dns_in)) {
            int64_t idx_in = ups_offset_in + dns_in_idx;
            auto [dns_flip, coeff] = term_action(dns_in);

            if constexpr (fermi_ups) { // not ideal to do this here
              if (bits::popcnt(ups_in) & 1) {
                coeff = -coeff;
              }
            }

            auto [idx_dns_flip, fermi_dn, sym] =
                basis_out.index_dns_fermi_sym(dns_flip, syms_out, dnss_out);

            coeff *= bloch_factors(sym);

            if (idx_dns_flip != invalid_index) {
              int64_t idx_out = ups_offset_out + idx_dns_flip;
              bool fermi_up = basis_out.fermi_bool_ups(sym, ups_out);
              coeff_t val =
                  coeff * norms_out[idx_dns_flip] / norms_in[dns_in_idx];

              XDIAG_FILL(idx_in, idx_out, (fermi_up ^ fermi_dn) ? -val : val);
            }
          }

          ++dns_in_idx;
        }
      } // if non trivial stabilizer
    } // loop over ups
#ifdef _OPENMP
  }
#endif
}

template <typename coeff_t, bool fermi_ups, class basis_t,
          class non_zero_term_f, class term_action_f, class fill_f>
void generic_term_dns_no_sym(basis_t const &basis_in, basis_t const &basis_out,
                             non_zero_term_f non_zero_term,
                             term_action_f term_action, fill_f fill) {
  using bit_t = typename basis_t::bit_t;

  int64_t size_ups_in = basis_in.size_ups();
  int64_t size_ups_out = basis_out.size_ups();

  int64_t size_dns_in = basis_in.size_dns();
  int64_t size_dns_out = basis_out.size_dns();
  assert(size_ups_in == size_ups_out);

#ifdef _OPENMP
#pragma omp parallel
  {
    int num_thread = omp_get_thread_num();
    auto dns_and_idces = basis_in.states_indices_dns_thread();
#else
  auto dns_and_idces = basis_in.states_indices_dns();
#endif
    for (auto [dns_in, idx_dns_in] : dns_and_idces) {
      if (non_zero_term(dns_in)) {
        auto [dns_out, coeff] = term_action(dns_in);

        int64_t idx_dns_out = basis_out.index_dns(dns_out);
        int64_t idx_out_start = idx_dns_out;
        int64_t idx_out_end = idx_dns_out + basis_out.size();

        if constexpr (fermi_ups) {
          int64_t idx_out = idx_out_start;
          int64_t idx_in = idx_dns_in;
          for (bit_t ups : basis_in.states_ups()) {
            bool fermi = bits::popcnt(ups) & 1;
            XDIAG_FILL(idx_in, idx_out, fermi ? coeff : -coeff);
            idx_out += size_dns_out;
            idx_in += size_dns_in;
          }
        } else {
          for (int64_t idx_out = idx_out_start, idx_in = idx_dns_in;
               idx_out < idx_out_end;
               idx_out += size_dns_out, idx_in += size_dns_in) {
            XDIAG_FILL(idx_in, idx_out, coeff);
          }
        }
      }
    }

#ifdef _OPENMP
  }
#endif
}
template <bool symmetric, typename coeff_t, bool fermi_ups, class basis_t,
          class non_zero_term_f, class term_action_f, class fill_f>
void generic_term_dns(basis_t const &basis_in, basis_t const &basis_out,
                      non_zero_term_f non_zero_term, term_action_f term_action,
                      fill_f fill) {
  if constexpr (symmetric) {
    generic_term_dns_sym(basis_in, basis_out, non_zero_term, term_action, fill);
  } else { // if not symmetric
    generic_term_dns_no_sym(basis_in, basis_out, non_zero_term, term_action,
                            fill);
  }
}

} // namespace xdiag::basis::electron
