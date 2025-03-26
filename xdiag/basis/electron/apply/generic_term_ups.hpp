#pragma once

#include <vector>

namespace xdiag::basis::electron {

template <typename bit_t, typename coeff_t, bool symmetric, class basis_t,
          class non_zero_term_f, class term_action_f, class fill_f>
void generic_term_ups(basis_t const &basis_in, basis_t const &basis_out,
                      non_zero_term_f non_zero_term, term_action_f term_action,
                      fill_f fill) {

  if constexpr (symmetric) {

    auto const &group_action = basis_out.group_action();
    Representation const &irrep = basis_out.irrep();
    auto bloch_factors = irrep.characters().as<arma::Col<coeff_t>>();

#ifdef _OPENMP
#pragma omp parallel for schedule(guided)
#endif
    for (int64_t idx_up_in = 0; idx_up_in < basis_in.n_rep_ups(); ++idx_up_in) {
      bit_t ups_in = basis_in.rep_ups(idx_up_in);
      if (non_zero_term(ups_in)) {

        auto [ups_flip, coeff] = term_action(ups_in);
        int64_t idx_ups_flip = basis_out.index_ups(ups_flip);
        bit_t ups_flip_rep = basis_out.rep_ups(idx_ups_flip);

        // Get limits, syms, and dns for ingoing ups
        int64_t ups_offset_in = basis_in.ups_offset(idx_up_in);
        auto dnss_in = basis_in.dns_for_ups_rep(ups_in);
        auto norms_in = basis_in.norms_for_ups_rep(ups_in);

        // Get limits, syms, and dns for outgoing ups
        int64_t ups_offset_out = basis_out.ups_offset(idx_ups_flip);
        auto syms_ups_out = basis_out.syms_ups(ups_flip);
        auto dnss_out = basis_out.dns_for_ups_rep(ups_flip_rep);
        auto norms_out = basis_out.norms_for_ups_rep(ups_flip_rep);

        // trivial up-stabilizer (likely)
        if (syms_ups_out.size() == 1) {
          int64_t sym = syms_ups_out.front();
          coeff_t prefac = coeff * bloch_factors(sym);
          bool fermi_up = basis_out.fermi_bool_ups(sym, ups_flip);

          int64_t idx_dn = 0;
          for (bit_t dns : dnss_in) {
            int64_t idx_in = ups_offset_in + idx_dn;
            bit_t dns_rep = group_action.apply(sym, dns);
            int64_t idx_out = ups_offset_out + basis_out.index_dns(dns_rep);
            bool fermi_dn = basis_out.fermi_bool_dns(sym, dns);

            coeff_t val =
                prefac / norms_in[idx_dn]; // norms_out = 1.0 in this case
            fill(idx_in, idx_out, (fermi_up ^ fermi_dn) ? -val : val);
            ++idx_dn;
          }

        } else { // non-trivial up-stabilizer (unlikely)
          auto syms = syms_ups_out;
          std::vector<coeff_t> prefacs(bloch_factors.size());
          for (int64_t i = 0; i < (int64_t)bloch_factors.size(); ++i) {
            prefacs[i] = coeff * bloch_factors(i);
          }

          int64_t idx_dn = 0;
          for (bit_t dns : dnss_in) {
            int64_t idx_in = ups_offset_in + idx_dn;

            auto [idx_dn_out, fermi_dn, sym] =
                basis_out.index_dns_fermi_sym(dns, syms, dnss_out);

            if (idx_dn_out != invalid_index) {
              int64_t idx_out = ups_offset_out + idx_dn_out;
              bool fermi_up = basis_out.fermi_bool_ups(sym, ups_flip);
              coeff_t val =
                  prefacs[sym] * norms_out[idx_dn_out] / norms_in[idx_dn];

              fill(idx_in, idx_out, (fermi_up ^ fermi_dn) ? -val : val);
            }
            ++idx_dn;
          }
        } // if trivial stabilizer or not
      } // if non_zero_term
    } // loop over ups

  } else { // not symmetric
    int64_t size_dns_in = basis_in.size_dns();
    int64_t size_dns_out = basis_out.size_dns();
    assert(size_dns_in == size_dns_out);

#ifdef _OPENMP
#pragma omp parallel
    {
      auto ups_and_idces = basis_in.states_indices_ups_thread();
#else
    auto ups_and_idces = basis_in.states_indices_ups();
#endif

      for (auto [ups_in, idx_ups_in] : ups_and_idces) {
        if (non_zero_term(ups_in)) {
          auto [ups_out, coeff] = term_action(ups_in);
          int64_t idx_ups_out = basis_out.index_ups(ups_out);

          int64_t idx_in_start = idx_ups_in * size_dns_in;
          int64_t idx_out_start = idx_ups_out * size_dns_out;
          int64_t idx_out_end = (idx_ups_out + 1) * size_dns_out;

          for (int64_t idx_out = idx_out_start, idx_in = idx_in_start;
               idx_out < idx_out_end; ++idx_out, ++idx_in) {
            fill(idx_in, idx_out, coeff);
          }
        }
      }

#ifdef _OPENMP
    }
#endif
  }
}

} // namespace xdiag::basis::electron
