#pragma once

#include <functional>

#include <xdiag/bits/bitops.hpp>

namespace xdiag::basis::tj {

template <typename coeff_t, bool symmetric, bool fermi_ups, class basis_t,
          class non_zero_term_ups_f, class non_zero_term_dns_f,
          class term_action_f, class fill_f>
void generic_term_dns(basis_t const &basis_in, basis_t const &basis_out,
                      non_zero_term_ups_f non_zero_term_ups,
                      non_zero_term_dns_f non_zero_term_dns,
                      term_action_f term_action, fill_f fill) {
  using bit_t = typename basis_t::bit_t;

  int64_t nsites = basis_in.nsites();
  assert(nsites == basis_out.nsites());
  bit_t sitesmask = ((bit_t)1 << nsites) - 1;

  if constexpr (symmetric) {

    Representation const &irrep = basis_out.irrep();
    auto bloch_factors = irrep.characters().as<arma::Col<coeff_t>>();

#ifdef _OPENMP
#pragma omp parallel for schedule(guided)
#endif
    for (int64_t idx_up_in = 0; idx_up_in < basis_in.n_rep_ups(); ++idx_up_in) {
      bit_t up_in = basis_in.rep_ups(idx_up_in);

      if (non_zero_term_ups(up_in)) {

        bit_t not_up_in = (~up_in) & sitesmask;
        int64_t up_in_offset = basis_in.ups_offset(idx_up_in);
        int64_t up_out_offset = basis_out.ups_offset(idx_up_in);

        auto dncs_in = basis_in.dns_for_ups_rep(up_in);
        auto syms_in = basis_in.syms_ups(up_in);
        auto norms_in = basis_in.norms_for_ups_rep(up_in);

        auto dncs_out = basis_out.dns_for_ups_rep(up_in);
        auto syms_out = basis_out.syms_ups(up_in);
        auto norms_out = basis_out.norms_for_ups_rep(up_in);

        // trivial stabilizer of up_in -> dns have to be deposited
        if (syms_in.size() == 1) {
          int64_t idx_in = up_in_offset;
          for (bit_t dnc_in : dncs_in) {
            int64_t dn_in = bits::deposit(dnc_in, not_up_in);
            if (non_zero_term_dns(dn_in)) {
              auto [dn_flip, coeff] = term_action(dn_in);
              if ((dn_flip & up_in) == 0) { // tJ constraint
                bit_t dnc_flip = bits::extract(dn_flip, not_up_in);
                int64_t idx_dnc_flip = basis_out.dnsc_index(dnc_flip);
                int64_t idx_out = up_out_offset + idx_dnc_flip;

                if constexpr (fermi_ups) {
                  bool fermi_up = (bool)(bits::popcnt(up_in) & 1);
                  fill(idx_in, idx_out, fermi_up ? -coeff : coeff);
                } else {
                  fill(idx_in, idx_out, coeff);
                }
              }
            }
            ++idx_in;
          }
        }

        // non-trivial stabilizer of ups -> dns don't have to be deposited
        else {
          int64_t idx_in = up_in_offset;
          int64_t idx_dn_in = 0;
          for (bit_t dn_in : dncs_in) {
            if (non_zero_term_dns(dn_in)) {
              auto [dn_flip, coeff] = term_action(dn_in);
              auto [idx_dn_flip, fermi_dn, sym] =
                  basis_out.index_dns_fermi_sym(dn_flip, syms_out, dncs_out);

              if (idx_dn_flip != invalid_index) {
                int64_t idx_out = up_out_offset + idx_dn_flip;
                bool fermi_up = basis_out.fermi_bool_ups(sym, up_in);
                if constexpr (fermi_ups) {
                  fermi_up ^= (bool)(bits::popcnt(up_in) & 1);
                }
                coeff_t val = coeff * bloch_factors(sym) *
                              norms_out[idx_dn_flip] / norms_in[idx_dn_in];
                fill(idx_in, idx_out, (fermi_up ^ fermi_dn) ? -val : val);
              }
            }
            ++idx_dn_in;
            ++idx_in;
          }
        }
      }
    }

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

          int64_t up_in_offset = basis_in.ups_offset(idx_up_in);
          int64_t up_out_offset = basis_out.ups_offset(idx_up_in);
          bit_t not_up_in = (~up_in) & sitesmask;

          auto dncs_in = basis_in.states_dncs(up_in);
          int64_t idx_in = up_in_offset;
          for (bit_t dnc_in : dncs_in) {
            bit_t dn_in = bits::deposit(dnc_in, not_up_in);

            if (non_zero_term_dns(dn_in)) {
              auto [dn_flip, coeff] = term_action(dn_in);
              if ((up_in & dn_flip) == 0) {
                bit_t dnc_flip = bits::extract(dn_flip, not_up_in);
                int64_t idx_dn_flip = basis_out.index_dncs(dnc_flip);
                int64_t idx_out = up_out_offset + idx_dn_flip;

                if constexpr (fermi_ups) {
                  bool fermi_up = (bool)(bits::popcnt(up_in) & 1);
                  fill(idx_in, idx_out, fermi_up ? -coeff : coeff);
                } else {
                  fill(idx_in, idx_out, coeff);
                }
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

} // namespace xdiag::basis::tj
