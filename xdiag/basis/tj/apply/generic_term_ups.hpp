// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <vector>

#include <xdiag/bits/bitops.hpp>
#include <xdiag/parallel/omp/omp_utils.hpp>

namespace xdiag::basis::tj {

template <typename coeff_t, bool symmetric, class basis_t,
          class non_zero_term_f, class term_action_f, class fill_f>
void generic_term_ups(basis_t const &basis_in, basis_t const &basis_out,
                      non_zero_term_f non_zero_term, term_action_f term_action,
                      fill_f fill) {
  using bit_t = typename basis_t::bit_t;

  int64_t nsites = basis_in.nsites();
  assert(nsites == basis_out.nsites());
  bit_t sitesmask = ((bit_t)1 << nsites) - 1;

  if constexpr (symmetric) {

    auto const &group_action = basis_out.group_action();
    Representation const &irrep = basis_out.irrep();
    auto bloch_factors = irrep.characters().as<arma::Col<coeff_t>>();

#ifdef _OPENMP
#pragma omp parallel
    {
      int num_thread = omp_get_thread_num();
#pragma omp for schedule(runtime)
#endif
      for (int64_t idx_up_in = 0; idx_up_in < basis_in.n_rep_ups();
           ++idx_up_in) {
        bit_t ups_in = basis_in.rep_ups(idx_up_in);
        if (non_zero_term(ups_in)) {

          auto [ups_flip, coeff] = term_action(ups_in);
          int64_t idx_ups_flip = basis_out.index_ups(ups_flip);
          bit_t ups_flip_rep = basis_out.rep_ups(idx_ups_flip);
          bit_t not_ups_flip_rep = (~ups_flip_rep) & sitesmask;

          // Get limits, syms, and dns for ingoing ups
          int64_t ups_offset_in = basis_in.ups_offset(idx_up_in);
          auto syms_ups_in = basis_in.syms_ups(ups_in);
          auto dnss_in = basis_in.dns_for_ups_rep(ups_in);
          auto norms_in = basis_in.norms_for_ups_rep(ups_in);

          // Get limits, syms, and dns for outgoing ups
          int64_t ups_offset_out = basis_out.ups_offset(idx_ups_flip);
          auto syms_ups_out = basis_out.syms_ups(ups_flip);
          auto dnss_out = basis_out.dns_for_ups_rep(ups_flip_rep);
          auto norms_out = basis_out.norms_for_ups_rep(ups_flip_rep);

          ////////////////////////////////////////////////////////////////////////
          // Trivial stabilizer of target ups
          if (syms_ups_out.size() == 1) {
            int64_t sym = syms_ups_out.front();
            coeff_t prefac = coeff * bloch_factors(sym);
            bool fermi_up = basis_out.fermi_bool_ups(sym, ups_flip);

            // Origin ups trivial stabilizer -> dns need to be deposited
            if (syms_ups_in.size() == 1) {
              int64_t idx_in = ups_offset_in;
              bit_t not_ups_in = (~ups_in) & sitesmask;
              for (bit_t dnsc : dnss_in) {
                bit_t dns = bits::deposit(dnsc, not_ups_in);
                if ((dns & ups_flip) == 0) { // t-J constraint
                  bit_t dns_rep = group_action.apply(sym, dns);
                  bit_t dns_rep_c = bits::extract(dns_rep, not_ups_flip_rep);
                  int64_t idx_out =
                      ups_offset_out + basis_out.dnsc_index(dns_rep_c);
                  bool fermi_dn = basis_out.fermi_bool_dns(sym, dns);
                  XDIAG_FILL(idx_in, idx_out,
                             (fermi_up ^ fermi_dn) ? -prefac : prefac);
                }
                ++idx_in;
              }
            }
            // Origin ups have stabilizer -> dns DONT need to be deposited
            else {
              int64_t idx_dn = 0;
              int64_t idx_in = ups_offset_in;
              for (bit_t dns : dnss_in) {
                if ((dns & ups_flip) == 0) { // t-J constraint
                  auto [idx_dn_out, fermi_dn] =
                      basis_out.index_dns_fermi(dns, sym, not_ups_flip_rep);
                  coeff_t val = prefac / norms_in[idx_dn];
                  int64_t idx_out = ups_offset_out + idx_dn_out;
                  XDIAG_FILL(idx_in, idx_out,
                             (fermi_up ^ fermi_dn) ? -val : val);
                }
                ++idx_in;
                ++idx_dn;
              }
            }

            ////////////////////////////////////////////////////////////////////////
            // Target ups have non-trivial stabilizer
          } else {
            std::vector<coeff_t> prefacs(bloch_factors.size());
            for (int64_t i = 0; i < (int64_t)bloch_factors.size(); ++i) {
              prefacs[i] = coeff * bloch_factors(i);
            }

            // Origin ups trivial stabilizer -> dns need to be deposited
            if (syms_ups_in.size() == 1) {
              bit_t not_ups_in = (~ups_in) & sitesmask;

              int64_t idx_in = ups_offset_in;
              for (bit_t dnsc : dnss_in) {
                bit_t dns = bits::deposit(dnsc, not_ups_in);

                if ((dns & ups_flip) == 0) { // t-J constraint
                  auto [idx_dn_out, fermi_dn, sym] =
                      basis_out.index_dns_fermi_sym(dns, syms_ups_out,
                                                    dnss_out);

                  if (idx_dn_out != invalid_index) {
                    int64_t idx_out = ups_offset_out + idx_dn_out;
                    bool fermi_up = basis_out.fermi_bool_ups(sym, ups_flip);
                    coeff_t val = prefacs[sym] * norms_out[idx_dn_out];

                    XDIAG_FILL(idx_in, idx_out,
                               (fermi_up ^ fermi_dn) ? -val : val);
                  }
                }
                ++idx_in;
              }
            }

            // Origin ups non-trivial stabilizer -> dns DONT need to be
            // deposited
            else {
              int64_t idx_in = ups_offset_in;
              int64_t idx_dn = 0;
              for (bit_t dns : dnss_in) {
                if ((dns & ups_flip) == 0) { // t-J constraint
                  auto [idx_dn_out, fermi_dn, sym] =
                      basis_out.index_dns_fermi_sym(dns, syms_ups_out,
                                                    dnss_out);
                  if (idx_dn_out != invalid_index) {
                    int64_t idx_out = ups_offset_out + idx_dn_out;
                    bool fermi_up = basis_out.fermi_bool_ups(sym, ups_flip);
                    coeff_t val =
                        prefacs[sym] * norms_out[idx_dn_out] / norms_in[idx_dn];
                    XDIAG_FILL(idx_in, idx_out,
                               (fermi_up ^ fermi_dn) ? -val : val);
                  }
                }
                ++idx_dn;
                ++idx_in;
              }
            }

          } // if target trivial stabilizer or not
        } // if non_zero_term
      } // loop over ups
#ifdef _OPENMP
    }
#endif

  } else { // if not symmetric
#ifdef _OPENMP
#pragma omp parallel
    {
      int num_thread = omp_get_thread_num();
      auto ups_and_idces = basis_in.states_indices_ups_thread();
#else
    auto ups_and_idces = basis_in.states_indices_ups();
#endif
      for (auto [up_in, idx_up_in] : ups_and_idces) {
        if (non_zero_term(up_in)) {

          auto [up_flip, coeff] = term_action(up_in);
          bit_t not_up_in = (~up_in) & sitesmask;
          bit_t not_up_flip = (~up_flip) & sitesmask;
          int64_t idx_up_flip = basis_out.index_ups(up_flip);
          int64_t idx_up_flip_offset = basis_out.ups_offset(idx_up_flip);

          auto dncs_in = basis_in.states_dncs(up_in);
          int64_t idx_in = basis_in.ups_offset(idx_up_in);
          for (bit_t dnc_in : dncs_in) {
            bit_t dn_in = bits::deposit(dnc_in, not_up_in);
            if ((up_flip & dn_in) == 0) { // tJ constraint
              bit_t dnc_out = bits::extract(dn_in, not_up_flip);
              int64_t idx_dnc_out = basis_out.index_dncs(dnc_out);
              int64_t idx_out = idx_up_flip_offset + idx_dnc_out;
              XDIAG_FILL(idx_in, idx_out, coeff);
            }
            ++idx_in;
          }
        }
      } // loop over ups
#ifdef _OPENMP
    }
#endif
  } // if not symmetric
}

} // namespace xdiag::basis::tj
