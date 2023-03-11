#pragma once

#include <vector>

namespace hydra::tj {

template <typename bit_t, typename coeff_t, bool symmetric, class IndexingIn,
          class IndexingOut, class NonZeroTermUps, class NonZeroTermDns,
          class TermActionUps, class TermActionDns, class Fill>
void generic_term_mixed(IndexingIn &&indexing_in, IndexingOut &&indexing_out,
                        NonZeroTermUps &&non_zero_term_ups,
                        NonZeroTermDns &&non_zero_term_dns,
                        TermActionUps &&term_action_ups,
                        TermActionDns &&term_action_dns, Fill &&fill) {
  int n_sites = indexing_in.n_sites();
  assert(n_sites == indexing_out.n_sites());
  bit_t sitesmask = ((bit_t)1 << n_sites) - 1;

  if constexpr (symmetric) {

    auto const &group_action = indexing_out.group_action();
    auto const &irrep = indexing_out.irrep();
    std::vector<coeff_t> bloch_factors;
    if constexpr (is_complex<coeff_t>()) {
      bloch_factors = irrep.characters();
    } else {
      bloch_factors = irrep.characters_real();
    }

    for (idx_t up_in_idx = 0; up_in_idx < indexing_in.n_rep_ups();
         ++up_in_idx) {
      bit_t up_in = indexing_in.rep_ups(up_in_idx);
      bit_t not_up_in = (~up_in) & sitesmask;

      if (non_zero_term_ups(up_in)) {

        auto [up_flip, coeff_up] = term_action_ups(up_in);
        idx_t idx_up_flip = indexing_out.index_ups(up_flip);
        bit_t up_flip_rep = indexing_out.rep_ups(idx_up_flip);
        bit_t not_up_flip_rep = (~up_flip_rep) & sitesmask;

        // Get limits, syms, and dns for ingoing ups
        idx_t up_in_offset = indexing_in.ups_offset(up_in_idx);
        auto up_in_syms = indexing_in.syms_ups(up_in);
        auto dnss_in = indexing_in.dns_for_ups_rep(up_in);
        auto norms_in = indexing_in.norms_for_ups_rep(up_in);

        // Get limits, syms, and dns for outgoing ups
        idx_t up_out_offset = indexing_out.ups_offset(idx_up_flip);
        auto up_out_syms = indexing_out.syms_ups(up_flip);
        auto dnss_out = indexing_out.dns_for_ups_rep(up_flip_rep);
        auto norms_out = indexing_out.norms_for_ups_rep(up_flip_rep);

        ////////////////////////////////////////////////////////////////////////
        // Trivial stabilizer of target ups
        if (up_out_syms.size() == 1) {
          int sym = up_out_syms.front();
          coeff_t prefac = coeff_up * bloch_factors[sym];
          bool fermi_up = indexing_out.fermi_bool_ups(sym, up_flip);

          // Origin ups trivial stabilizer -> dns need to be deposited
          if (up_in_syms.size() == 1) {
            idx_t idx_in = up_in_offset;
            for (bit_t dnc_in : dnss_in) {
              bit_t dn_in = bitops::deposit(dnc_in, not_up_in);
              if (non_zero_term_dns(dn_in)) {
                auto [dn_flip, coeff_dn] = term_action_dns(dn_in);

                if ((dn_flip & up_flip) == 0) { // t-J constraint
                  bit_t dn_rep = group_action.apply(sym, dn_flip);
                  bit_t dnc_rep = bitops::extract(dn_rep, not_up_flip_rep);
                  idx_t dnc_rep_idx = indexing_out.dnsc_index(dnc_rep);
                  idx_t idx_out = up_out_offset + dnc_rep_idx;
                  bool fermi_dn = indexing_out.fermi_bool_dns(sym, dn_flip);
                  fill(idx_out, idx_in,
                       (fermi_up ^ fermi_dn) ? -prefac * coeff_dn
                                             : prefac * coeff_dn);
                }
              }
              ++idx_in;
            }
          }
          // Origin ups have stabilizer -> dns DONT need to be deposited
          else {
            idx_t idx_in = up_in_offset;
            idx_t idx_dn = 0;
            for (bit_t dn : dnss_in) {
              if (non_zero_term_dns(dn)) {
                auto [dn_flip, coeff_dn] = term_action_dns(dn);
                if ((dn_flip & up_flip) == 0) { // t-J constraint
                  auto [idx_dn_out, fermi_dn] = indexing_out.index_dns_fermi(
                      dn_flip, sym, not_up_flip_rep);
                  coeff_t val = prefac * coeff_dn / norms_in[idx_dn];
                  idx_t idx_out = up_out_offset + idx_dn_out;
                  fill(idx_out, idx_in, (fermi_up ^ fermi_dn) ? -val : val);
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
          for (int i = 0; i < (int)bloch_factors.size(); ++i) {
            prefacs[i] = coeff_up * bloch_factors[i];
          }

          // Origin ups trivial stabilizer -> dns need to be deposited
          if (up_in_syms.size() == 1) {
            idx_t idx_in = up_in_offset;
            for (bit_t dnc : dnss_in) {
              bit_t dn = bitops::deposit(dnc, not_up_in);
              if (non_zero_term_dns(dn)) {
                auto [dn_flip, coeff_dn] = term_action_dns(dn);

                if ((dn_flip & up_flip) == 0) { // t-J constraint
                  auto [idx_dn_out, fermi_dn, sym] =
                      indexing_out.index_dns_fermi_sym(dn_flip, up_out_syms,
                                                       dnss_out);

                  if (idx_dn_out != invalid_index) {
                    idx_t idx_out = up_out_offset + idx_dn_out;
                    bool fermi_up = indexing_out.fermi_bool_ups(sym, up_flip);
                    coeff_t val =
                        prefacs[sym] * coeff_dn * norms_out[idx_dn_out];
                    fill(idx_out, idx_in, (fermi_up ^ fermi_dn) ? -val : val);
                  }
                }
              }
              ++idx_in;
            }
          }

          // Origin ups non-trivial stabilizer -> dns DONT need to be
          // deposited
          else {
            idx_t idx_in = up_in_offset;
            idx_t idx_dn = 0;
            for (bit_t dn : dnss_in) {
              if (non_zero_term_dns(dn)) {
                auto [dn_flip, coeff_dn] = term_action_dns(dn);
                if ((dn_flip & up_flip) == 0) { // t-J constraint
                  auto [idx_dn_out, fermi_dn, sym] =
                      indexing_out.index_dns_fermi_sym(dn_flip, up_out_syms,
                                                       dnss_out);
                  if (idx_dn_out != invalid_index) {
                    idx_t idx_out = up_out_offset + idx_dn_out;
                    bool fermi_up = indexing_out.fermi_bool_ups(sym, up_flip);
                    coeff_t val =
                        prefacs[sym] * norms_out[idx_dn_out] / norms_in[idx_dn];
                    fill(idx_out, idx_in, (fermi_up ^ fermi_dn) ? -val : val);
                  }
                }
                ++idx_in;
                ++idx_dn;
              }
            }
          }

        } // if target trivial stabilizer or not
      }   // if non_zero_term_ups
    }     // loop over ups

  } else { // if not symmetric

    idx_t size_dncs_in = indexing_in.size_dncs(); // wouldn't be valid for NoNp
    idx_t size_dncs_out = indexing_out.size_dncs();
    assert(size_dncs_in == size_dncs_out);

    auto ups_and_idces = indexing_in.states_indices_ups();
    for (auto [up_in, idx_up_in] : ups_and_idces) {
      if (non_zero_term_ups(up_in)) {

        auto [up_flip, coeff_up] = term_action_ups(up_in);
        bit_t not_up_in = (~up_in) & sitesmask;
        bit_t not_up_flip = (~up_flip) & sitesmask;
        idx_t idx_up_flip = indexing_out.index_ups(up_flip);
        idx_t idx_up_flip_offset = idx_up_flip * size_dncs_out;

        auto dncs_in = indexing_in.states_dncs(up_in);
        idx_t idx_in =
            idx_up_in * size_dncs_in; // wouldn't be valid for NoNp (use offset)
        for (bit_t dnc_in : dncs_in) {
          bit_t dn_in = bitops::deposit(dnc_in, not_up_in);
          if (non_zero_term_dns(dn_in)) {
            auto [dn_flip, coeff_dn] = term_action_dns(dn_in);
            if ((dn_flip & up_flip) == 0) { // t-J constraint
              bit_t dnc_out = bitops::extract(dn_in, not_up_flip);
              idx_t idx_dnc_out = indexing_out.index_dncs(dnc_out);
              idx_t idx_out = idx_up_flip_offset + idx_dnc_out;
              fill(idx_out, idx_in, coeff_up * coeff_dn);
            }
          }
          ++idx_in;
        }
      }
    } // loop over ups
  }   // if not symmetric
}

} // namespace hydra::tj
