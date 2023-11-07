#pragma once

#include <hydra/bits/bitops.h>
#include <vector>

namespace hydra::tj {

template <typename bit_t, typename coeff_t, bool symmetric, class BasisIn,
          class BasisOut, class NonZeroTermUps, class NonZeroTermDns,
          class TermActionUps, class TermActionDns, class Fill>
void generic_term_mixed(BasisIn &&basis_in, BasisOut &&basis_out,
                        NonZeroTermUps &&non_zero_term_ups,
                        NonZeroTermDns &&non_zero_term_dns,
                        TermActionUps &&term_action_ups,
                        TermActionDns &&term_action_dns, Fill &&fill) {
  int64_t n_sites = basis_in.n_sites();
  assert(n_sites == basis_out.n_sites());
  bit_t sitesmask = ((bit_t)1 << n_sites) - 1;

  if constexpr (symmetric) {

    // auto const &group_action = basis_out.group_action();
    auto const &irrep = basis_out.irrep();
    std::vector<coeff_t> bloch_factors;
    if constexpr (iscomplex<coeff_t>()) {
      bloch_factors = irrep.characters();
    } else {
      bloch_factors = irrep.characters_real();
    }

#ifdef _OPENMP
#pragma omp parallel for schedule(guided)
#endif
    for (int64_t up_in_idx = 0; up_in_idx < basis_in.n_rep_ups(); ++up_in_idx) {
      bit_t up_in = basis_in.rep_ups(up_in_idx);
      bit_t not_up_in = (~up_in) & sitesmask;

      if (non_zero_term_ups(up_in)) {

        auto [up_flip, coeff_up] = term_action_ups(up_in);
        int64_t idx_up_flip = basis_out.index_ups(up_flip);
        bit_t up_flip_rep = basis_out.rep_ups(idx_up_flip);
        bit_t not_up_flip_rep = (~up_flip_rep) & sitesmask;

        // Get limits, syms, and dns for ingoing ups
        int64_t up_in_offset = basis_in.ups_offset(up_in_idx);
        auto up_in_syms = basis_in.syms_ups(up_in);
        auto dnss_in = basis_in.dns_for_ups_rep(up_in);
        auto norms_in = basis_in.norms_for_ups_rep(up_in);

        // Get limits, syms, and dns for outgoing ups
        int64_t up_out_offset = basis_out.ups_offset(idx_up_flip);
        auto up_out_syms = basis_out.syms_ups(up_flip);
        auto dnss_out = basis_out.dns_for_ups_rep(up_flip_rep);
        auto norms_out = basis_out.norms_for_ups_rep(up_flip_rep);

        ////////////////////////////////////////////////////////////////////////
        // Trivial stabilizer of target ups
        if (up_out_syms.size() == 1) {
          int64_t sym = up_out_syms.front();
          coeff_t prefac = coeff_up * bloch_factors[sym];
          bool fermi_up = basis_out.fermi_bool_ups(sym, up_flip);

          // Origin ups trivial stabilizer -> dns need to be deposited
          if (up_in_syms.size() == 1) {
            int64_t idx_in = up_in_offset;
            for (bit_t dnc_in : dnss_in) {
              bit_t dn_in = bits::deposit(dnc_in, not_up_in);
              if (non_zero_term_dns(dn_in)) {
                auto [dn_flip, coeff_dn] = term_action_dns(dn_in);

                if ((dn_flip & up_flip) == 0) { // t-J constraint
                  auto [idx_dn_out, fermi_dn] =
                      basis_out.index_dns_fermi(dn_flip, sym, not_up_flip_rep);
                  int64_t idx_out = up_out_offset + idx_dn_out;
                  fill(idx_in, idx_out,
                       (fermi_up ^ fermi_dn) ? -prefac * coeff_dn
                                             : prefac * coeff_dn);
                }
              }
              ++idx_in;
            }
          }

          // Origin ups have stabilizer -> dns DONT need to be deposited
          else {
            int64_t idx_in = up_in_offset;
            int64_t idx_dn = 0;
            for (bit_t dn : dnss_in) {
              if (non_zero_term_dns(dn)) {
                auto [dn_flip, coeff_dn] = term_action_dns(dn);
                if ((dn_flip & up_flip) == 0) { // t-J constraint
                  auto [idx_dn_out, fermi_dn] =
                      basis_out.index_dns_fermi(dn_flip, sym, not_up_flip_rep);
                  coeff_t val = prefac * coeff_dn / norms_in[idx_dn];
                  int64_t idx_out = up_out_offset + idx_dn_out;
                  fill(idx_in, idx_out, (fermi_up ^ fermi_dn) ? -val : val);
                }
              }
              ++idx_in;
              ++idx_dn;
            }
          }
        }

        ////////////////////////////////////////////////////////////////////////
        // Target ups have non-trivial stabilizer
        else {
          std::vector<coeff_t> prefacs(bloch_factors.size());
          for (int64_t i = 0; i < (int64_t)bloch_factors.size(); ++i) {
            prefacs[i] = coeff_up * bloch_factors[i];
          }

          // Origin ups trivial stabilizer -> dns need to be deposited
          if (up_in_syms.size() == 1) {
            int64_t idx_in = up_in_offset;
            for (bit_t dnc : dnss_in) {
              bit_t dn = bits::deposit(dnc, not_up_in);
              if (non_zero_term_dns(dn)) {
                auto [dn_flip, coeff_dn] = term_action_dns(dn);

                if ((dn_flip & up_flip) == 0) { // t-J constraint
                  auto [idx_dn_out, fermi_dn, sym] =
                      basis_out.index_dns_fermi_sym(dn_flip, up_out_syms,
                                                    dnss_out);
                  if (idx_dn_out != invalid_index) {
                    int64_t idx_out = up_out_offset + idx_dn_out;
                    bool fermi_up = basis_out.fermi_bool_ups(sym, up_flip);
                    coeff_t val =
                        prefacs[sym] * coeff_dn * norms_out[idx_dn_out];
                    fill(idx_in, idx_out, (fermi_up ^ fermi_dn) ? -val : val);
                  }
                }
              }
              ++idx_in;
            }
          }

          // Origin ups non-trivial stabilizer -> dns DONT need to be
          // deposited
          else {

            int64_t idx_in = up_in_offset;
            int64_t idx_dn = 0;
            for (bit_t dn : dnss_in) {
              if (non_zero_term_dns(dn)) {
                auto [dn_flip, coeff_dn] = term_action_dns(dn);

                if ((dn_flip & up_flip) == 0) { // t-J constraint
                  auto [idx_dn_out, fermi_dn, sym] =
                      basis_out.index_dns_fermi_sym(dn_flip, up_out_syms,
                                                    dnss_out);
                  if (idx_dn_out != invalid_index) {
                    int64_t idx_out = up_out_offset + idx_dn_out;
                    bool fermi_up = basis_out.fermi_bool_ups(sym, up_flip);
                    coeff_t val = prefacs[sym] * coeff_dn *
                                  norms_out[idx_dn_out] / norms_in[idx_dn];
                    fill(idx_in, idx_out, (fermi_up ^ fermi_dn) ? -val : val);
                  }
                }
              }
              ++idx_in;
              ++idx_dn;
            }
          }

        } // if target trivial stabilizer or not
      }   // if non_zero_term_ups
    }     // loop over ups

  } else { // if not symmetric

#ifdef _OPENMP
#pragma omp parallel
    {
      auto ups_and_idces = basis_in.states_indices_ups_thread();
#else
    auto ups_and_idces = basis_in.states_indices_ups();
#endif
      for (auto [up_in, idx_up_in] : ups_and_idces) {
        if (non_zero_term_ups(up_in)) {

          auto [up_flip, coeff_up] = term_action_ups(up_in);
          bit_t not_up_in = (~up_in) & sitesmask;
          bit_t not_up_flip = (~up_flip) & sitesmask;
          int64_t idx_up_flip = basis_out.index_ups(up_flip);
          int64_t idx_up_flip_offset = basis_out.ups_offset(idx_up_flip);

          auto dncs_in = basis_in.states_dncs(up_in);
          int64_t idx_in = basis_in.ups_offset(idx_up_in);

          for (bit_t dnc_in : dncs_in) {
            bit_t dn_in = bits::deposit(dnc_in, not_up_in);
            if (non_zero_term_dns(dn_in)) {
              auto [dn_flip, coeff_dn] = term_action_dns(dn_in);
              if ((dn_flip & up_flip) == 0) { // t-J constraint
                bit_t dnc_out = bits::extract(dn_flip, not_up_flip);
                int64_t idx_dnc_out = basis_out.index_dncs(dnc_out);
                int64_t idx_out = idx_up_flip_offset + idx_dnc_out;
                fill(idx_in, idx_out, coeff_up * coeff_dn);
                // Log("tJ in: {} out: {} coeff: {}", idx_in, idx_out,
                //     coeff_up * coeff_dn);
              }
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

} // namespace hydra::tj
