#pragma once

#include <vector>
#include <hydra/bits/bitops.h>

namespace hydra::tj {

template <typename bit_t, typename coeff_t, bool symmetric, class BasisIn,
          class BasisOut, class NonZeroTerm, class TermAction, class Fill>
void generic_term_ups(BasisIn &&basis_in, BasisOut &&basis_out,
                      NonZeroTerm &&non_zero_term, TermAction &&term_action,
                      Fill &&fill) {
  int64_t n_sites = basis_in.n_sites();
  assert(n_sites == basis_out.n_sites());
  bit_t sitesmask = ((bit_t)1 << n_sites) - 1;

  if constexpr (symmetric) {

    auto const &group_action = basis_out.group_action();
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
    for (idx_t idx_up_in = 0; idx_up_in < basis_in.n_rep_ups(); ++idx_up_in) {
      bit_t ups_in = basis_in.rep_ups(idx_up_in);
      if (non_zero_term(ups_in)) {

        auto [ups_flip, coeff] = term_action(ups_in);
        idx_t idx_ups_flip = basis_out.index_ups(ups_flip);
        bit_t ups_flip_rep = basis_out.rep_ups(idx_ups_flip);
        bit_t not_ups_flip_rep = (~ups_flip_rep) & sitesmask;

        // Get limits, syms, and dns for ingoing ups
        idx_t ups_offset_in = basis_in.ups_offset(idx_up_in);
        auto syms_ups_in = basis_in.syms_ups(ups_in);
        auto dnss_in = basis_in.dns_for_ups_rep(ups_in);
        auto norms_in = basis_in.norms_for_ups_rep(ups_in);

        // Get limits, syms, and dns for outgoing ups
        idx_t ups_offset_out = basis_out.ups_offset(idx_ups_flip);
        auto syms_ups_out = basis_out.syms_ups(ups_flip);
        auto dnss_out = basis_out.dns_for_ups_rep(ups_flip_rep);
        auto norms_out = basis_out.norms_for_ups_rep(ups_flip_rep);

        ////////////////////////////////////////////////////////////////////////
        // Trivial stabilizer of target ups
        if (syms_ups_out.size() == 1) {
          int64_t sym = syms_ups_out.front();
          coeff_t prefac = coeff * bloch_factors[sym];
          bool fermi_up = basis_out.fermi_bool_ups(sym, ups_flip);

          // Origin ups trivial stabilizer -> dns need to be deposited
          if (syms_ups_in.size() == 1) {
            idx_t idx_in = ups_offset_in;
            bit_t not_ups_in = (~ups_in) & sitesmask;
            for (bit_t dnsc : dnss_in) {
              bit_t dns = bits::deposit(dnsc, not_ups_in);
              if ((dns & ups_flip) == 0) { // t-J constraint
                bit_t dns_rep = group_action.apply(sym, dns);
                bit_t dns_rep_c = bits::extract(dns_rep, not_ups_flip_rep);
                idx_t idx_out =
                    ups_offset_out + basis_out.dnsc_index(dns_rep_c);
                bool fermi_dn = basis_out.fermi_bool_dns(sym, dns);
                fill(idx_out, idx_in, (fermi_up ^ fermi_dn) ? -prefac : prefac);
              }
              ++idx_in;
            }
          }
          // Origin ups have stabilizer -> dns DONT need to be deposited
          else {
            idx_t idx_dn = 0;
            idx_t idx_in = ups_offset_in;
            for (bit_t dns : dnss_in) {
              if ((dns & ups_flip) == 0) { // t-J constraint
                auto [idx_dn_out, fermi_dn] =
                    basis_out.index_dns_fermi(dns, sym, not_ups_flip_rep);
                coeff_t val = prefac / norms_in[idx_dn];
                idx_t idx_out = ups_offset_out + idx_dn_out;
                fill(idx_out, idx_in, (fermi_up ^ fermi_dn) ? -val : val);
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
            prefacs[i] = coeff * bloch_factors[i];
          }

          // Origin ups trivial stabilizer -> dns need to be deposited
          if (syms_ups_in.size() == 1) {
            bit_t not_ups_in = (~ups_in) & sitesmask;

            idx_t idx_in = ups_offset_in;
            for (bit_t dnsc : dnss_in) {
              bit_t dns = bits::deposit(dnsc, not_ups_in);

              if ((dns & ups_flip) == 0) { // t-J constraint
                auto [idx_dn_out, fermi_dn, sym] =
                    basis_out.index_dns_fermi_sym(dns, syms_ups_out, dnss_out);

                if (idx_dn_out != invalid_index) {
                  idx_t idx_out = ups_offset_out + idx_dn_out;
                  bool fermi_up = basis_out.fermi_bool_ups(sym, ups_flip);
                  coeff_t val = prefacs[sym] * norms_out[idx_dn_out];

                  fill(idx_out, idx_in, (fermi_up ^ fermi_dn) ? -val : val);
                }
              }
              ++idx_in;
            }
          }

          // Origin ups non-trivial stabilizer -> dns DONT need to be deposited
          else {
            idx_t idx_in = ups_offset_in;
            idx_t idx_dn = 0;
            for (bit_t dns : dnss_in) {
              if ((dns & ups_flip) == 0) { // t-J constraint
                auto [idx_dn_out, fermi_dn, sym] =
                    basis_out.index_dns_fermi_sym(dns, syms_ups_out, dnss_out);
                if (idx_dn_out != invalid_index) {
                  idx_t idx_out = ups_offset_out + idx_dn_out;
                  bool fermi_up = basis_out.fermi_bool_ups(sym, ups_flip);
                  coeff_t val =
                      prefacs[sym] * norms_out[idx_dn_out] / norms_in[idx_dn];
                  fill(idx_out, idx_in, (fermi_up ^ fermi_dn) ? -val : val);
                }
              }
              ++idx_dn;
              ++idx_in;
            }
          }

        } // if target trivial stabilizer or not
      }   // if non_zero_term
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
        if (non_zero_term(up_in)) {

          auto [up_flip, coeff] = term_action(up_in);
          bit_t not_up_in = (~up_in) & sitesmask;
          bit_t not_up_flip = (~up_flip) & sitesmask;
          idx_t idx_up_flip = basis_out.index_ups(up_flip);
          idx_t idx_up_flip_offset = basis_out.ups_offset(idx_up_flip);

          auto dncs_in = basis_in.states_dncs(up_in);
          idx_t idx_in = basis_in.ups_offset(idx_up_in);
          for (bit_t dnc_in : dncs_in) {
            bit_t dn_in = bits::deposit(dnc_in, not_up_in);
            if ((up_flip & dn_in) == 0) { // tJ constraint
              bit_t dnc_out = bits::extract(dn_in, not_up_flip);
              idx_t idx_dnc_out = basis_out.index_dncs(dnc_out);
              idx_t idx_out = idx_up_flip_offset + idx_dnc_out;
              fill(idx_out, idx_in, coeff);
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
