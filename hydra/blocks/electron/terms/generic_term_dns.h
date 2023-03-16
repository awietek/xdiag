#pragma once

#include <functional>

namespace hydra::electron {

template <typename bit_t, typename coeff_t, bool symmetric, bool fermi_ups,
          class IndexingIn, class IndexingOut, class NonZeroTerm,
          class TermAction, class Fill>
void generic_term_dns(IndexingIn &&indexing_in, IndexingOut &&indexing_out,
                      NonZeroTerm &&non_zero_term, TermAction &&term_action,
                      Fill &&fill) {

  if constexpr (symmetric) {

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
    for (idx_t idx_up = 0; idx_up < indexing_in.n_rep_ups(); ++idx_up) {
      bit_t ups_in = indexing_in.rep_ups(idx_up);
      bit_t ups_out = indexing_out.rep_ups(idx_up);
      assert(ups_out == ups_in);
      bit_t ups = ups_in;

      idx_t ups_offset_in = indexing_out.ups_offset(idx_up);
      auto dnss_in = indexing_in.dns_for_ups_rep(ups_in);
      auto norms_in = indexing_in.norms_for_ups_rep(ups_in);

      idx_t ups_offset_out = indexing_out.ups_offset(idx_up);
      auto syms_out = indexing_out.syms_ups(ups_out);
      auto dnss_out = indexing_out.dns_for_ups_rep(ups_out);
      auto norms_out = indexing_out.norms_for_ups_rep(ups_out);

      // trivial up-stabilizer (likely)
      if (syms_out.size() == 1) {
        idx_t dns_in_idx = 0;
        for (bit_t dns_in : dnss_in) {

          if (non_zero_term(dns_in)) {
            idx_t idx_in = ups_offset_in + dns_in_idx;
            auto [dns_flip, coeff] = term_action(dns_in);

            if constexpr (fermi_ups) { // not ideal to do this here
              if (bitops::popcnt(ups) & 1) {
                coeff = -coeff;
              }
            }

            idx_t idx_dns_flip = indexing_out.index_dns(dns_flip);
            idx_t idx_out = ups_offset_out + idx_dns_flip;
            fill(idx_out, idx_in, coeff);
          }
          ++dns_in_idx;
        }
      }

      else { // non-trivial up-stabilizer (unlikely)

        idx_t dns_in_idx = 0;
        for (bit_t dns_in : dnss_in) {

          if (non_zero_term(dns_in)) {
            idx_t idx_in = ups_offset_in + dns_in_idx;
            auto [dns_flip, coeff] = term_action(dns_in);

            if constexpr (fermi_ups) { // not ideal to do this here
              if (bitops::popcnt(ups) & 1) {
                coeff = -coeff;
              }
            }

            auto [idx_dns_flip, fermi_dn, sym] =
                indexing_out.index_dns_fermi_sym(dns_flip, syms_out, dnss_out);

            coeff *= bloch_factors[sym];

            if (idx_dns_flip != invalid_index) {
              idx_t idx_out = ups_offset_out + idx_dns_flip;
              bool fermi_up = indexing_out.fermi_bool_ups(sym, ups_out);
              coeff_t val =
                  coeff * norms_out[idx_dns_flip] / norms_in[dns_in_idx];

              fill(idx_out, idx_in, (fermi_up ^ fermi_dn) ? -val : val);
            }
          }

          ++dns_in_idx;
        }
      } // if non trivial stabilizer
    }   // loop over ups
  }

  else { // if not symmetric

    idx_t size_ups_in = indexing_in.size_ups();
    idx_t size_ups_out = indexing_out.size_ups();

    idx_t size_dns_in = indexing_in.size_dns();
    idx_t size_dns_out = indexing_out.size_dns();
    if (size_ups_in != size_ups_out) {
      Log.err("Error in apply_term_dns_no_sym: size of upspins space different "
              "for input and output indexing");
    }

#ifdef _OPENMP
#pragma omp parallel
    {
      auto dns_and_idces = indexing_in.states_indices_dns_thread();
#else
    auto dns_and_idces = indexing_in.states_indices_dns();
#endif
      for (auto [dns_in, idx_dns_in] : dns_and_idces) {
        if (non_zero_term(dns_in)) {
          auto [dns_out, coeff] = term_action(dns_in);

          idx_t idx_dns_out = indexing_out.index_dns(dns_out);
          idx_t idx_out_start = idx_dns_out;
          idx_t idx_out_end = idx_dns_out + indexing_out.size();

          if constexpr (fermi_ups) {
            idx_t idx_out = idx_out_start;
            idx_t idx_in = idx_dns_in;
            for (bit_t ups : indexing_in.states_ups()) {
              bool fermi = bitops::popcnt(ups) & 1;
              fill(idx_out, idx_in, fermi ? coeff : -coeff);
              idx_out += size_dns_out;
              idx_in += size_dns_in;
            }
          } else {
            for (idx_t idx_out = idx_out_start, idx_in = idx_dns_in;
                 idx_out < idx_out_end;
                 idx_out += size_dns_out, idx_in += size_dns_in) {
              fill(idx_out, idx_in, coeff);
            }
          }
        }
      }

#ifdef _OPENMP
    }
#endif
  }
}

} // namespace hydra::electron
