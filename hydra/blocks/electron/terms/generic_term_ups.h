#pragma once

#include <vector>

namespace hydra::electron {

template <typename bit_t, typename coeff_t, bool symmetric, class IndexingIn,
          class IndexingOut, class NonZeroTerm, class TermAction, class Fill>
void generic_term_ups(IndexingIn &&indexing_in, IndexingOut &&indexing_out,
                      NonZeroTerm &&non_zero_term, TermAction &&term_action,
                      Fill &&fill) {

  if constexpr (symmetric) {

    auto const &group_action = indexing_out.group_action();
    auto const &irrep = indexing_out.irrep();
    std::vector<coeff_t> bloch_factors;
    if constexpr (is_complex<coeff_t>()) {
      bloch_factors = irrep.characters();
    } else {
      bloch_factors = irrep.characters_real();
    }

#ifdef _OPENMP
#pragma omp parallel for schedule(guided)
#endif
    for (idx_t idx_up_in = 0; idx_up_in < indexing_in.n_rep_ups();
         ++idx_up_in) {
      bit_t ups_in = indexing_in.rep_ups(idx_up_in);
      if (non_zero_term(ups_in)) {

        auto [ups_flip, coeff] = term_action(ups_in);
        idx_t idx_ups_flip = indexing_out.index_ups(ups_flip);
        bit_t ups_flip_rep = indexing_out.rep_ups(idx_ups_flip);

        // Get limits, syms, and dns for ingoing ups
        idx_t ups_offset_in = indexing_in.ups_offset(idx_up_in);
        auto dnss_in = indexing_in.dns_for_ups_rep(ups_in);
        auto norms_in = indexing_in.norms_for_ups_rep(ups_in);

        // Get limits, syms, and dns for outgoing ups
        idx_t ups_offset_out = indexing_out.ups_offset(idx_ups_flip);
        auto syms_ups_out = indexing_out.syms_ups(ups_flip);
        auto dnss_out = indexing_out.dns_for_ups_rep(ups_flip_rep);
        auto norms_out = indexing_out.norms_for_ups_rep(ups_flip_rep);

        // trivial up-stabilizer (likely)
        if (syms_ups_out.size() == 1) {
          int sym = syms_ups_out.front();
          coeff_t prefac = coeff * bloch_factors[sym];
          bool fermi_up = indexing_out.fermi_bool_ups(sym, ups_flip);

          idx_t idx_dn = 0;
          for (bit_t dns : dnss_in) {
            idx_t idx_in = ups_offset_in + idx_dn;
            bit_t dns_rep = group_action.apply(sym, dns);
            idx_t idx_out = ups_offset_out + indexing_out.index_dns(dns_rep);
            bool fermi_dn = indexing_out.fermi_bool_dns(sym, dns);

            coeff_t val =
                prefac / norms_in[idx_dn]; // norms_out = 1.0 in this case
            fill(idx_out, idx_in, (fermi_up ^ fermi_dn) ? -val : val);
            ++idx_dn;
          }

        } else { // non-trivial up-stabilizer (unlikely)
          auto syms = syms_ups_out;
          std::vector<coeff_t> prefacs(bloch_factors.size());
          for (int i = 0; i < (int)bloch_factors.size(); ++i) {
            prefacs[i] = coeff * bloch_factors[i];
          }

          idx_t idx_dn = 0;
          for (bit_t dns : dnss_in) {
            idx_t idx_in = ups_offset_in + idx_dn;

            auto [idx_dn_out, fermi_dn, sym] =
                indexing_out.index_dns_fermi_sym(dns, syms, dnss_out);

            if (idx_dn_out != invalid_index) {
              idx_t idx_out = ups_offset_out + idx_dn_out;
              bool fermi_up = indexing_out.fermi_bool_ups(sym, ups_flip);
              coeff_t val =
                  prefacs[sym] * norms_out[idx_dn_out] / norms_in[idx_dn];

              fill(idx_out, idx_in, (fermi_up ^ fermi_dn) ? -val : val);
            }
            ++idx_dn;
          }
        } // if trivial stabilizer or not
      }   // if non_zero_term
    }     // loop over ups

  } else { // not symmetric
    idx_t size_dns_in = indexing_in.size_dns();
    idx_t size_dns_out = indexing_out.size_dns();
    if (size_dns_in != size_dns_out) {
      Log.err("Error in apply_term_ups: size of dnspins space different "
              "for input and output indexing");
    }

#ifdef _OPENMP
#pragma omp parallel
    {
      auto ups_and_idces = indexing_in.states_indices_ups_thread();
#else
    auto ups_and_idces = indexing_in.states_indices_ups();
#endif

      for (auto [ups_in, idx_ups_in] : ups_and_idces) {
        if (non_zero_term(ups_in)) {
          auto [ups_out, coeff] = term_action(ups_in);
          idx_t idx_ups_out = indexing_out.index_ups(ups_out);

          idx_t idx_in_start = idx_ups_in * size_dns_in;
          idx_t idx_out_start = idx_ups_out * size_dns_out;
          idx_t idx_out_end = (idx_ups_out + 1) * size_dns_out;

          for (idx_t idx_out = idx_out_start, idx_in = idx_in_start;
               idx_out < idx_out_end; ++idx_out, ++idx_in) {
            fill(idx_out, idx_in, coeff);
          }
        }
      }

#ifdef _OPENMP
    }
#endif
  }
}

} // namespace hydra::electron
