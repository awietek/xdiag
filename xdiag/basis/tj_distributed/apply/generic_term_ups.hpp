#pragma once

#include <vector>

#include <xdiag/bits/bitops.hpp>
#include <xdiag/utils/logger.hpp>

namespace xdiag::basis::tj_distributed {

template <typename bit_t, typename coeff_t, class basis_t,
          class non_zero_term_ups_f, class non_zero_term_dns_f,
          class term_action_f>
void generic_term_ups(basis_t const &basis_in, basis_t const &basis_out,
                      non_zero_term_ups_f non_zero_term_ups,
                      non_zero_term_dns_f non_zero_term_dns,
                      term_action_f term_action, const coeff_t *vec_in,
                      coeff_t *vec_out) {
  int64_t nsites = basis_in.nsites();
  assert(nsites == basis_out.nsites());
  bit_t sitesmask = ((bit_t)1 << nsites) - 1;

  int64_t nup_in = basis_in.nup();
  int64_t ndn_in = basis_in.ndn();
  int64_t nup_configurations_in =
      combinatorics::binomial(nsites - ndn_in, nup_in);

  int64_t nup_out = basis_out.nup();
  int64_t ndn_out = basis_out.ndn();
  int64_t nup_configurations_out =
      combinatorics::binomial(nsites - ndn_out, nup_out);

  // Loop over all configurations
  int64_t idx_dn = 0;
  for (bit_t dn : basis_in.my_dns()) {

    if (non_zero_term_dns(dn)) {
      bit_t not_dn = (~dn) & sitesmask;
      int64_t dn_offset_in = idx_dn * nup_configurations_in;

      for (int64_t idx_in = dn_offset_in;
           idx_in < dn_offset_in + nup_configurations_in; ++idx_in) {

        bit_t up = basis_in.my_ups_for_dns_storage(idx_in);

        // Check if hopping is possible
        if (non_zero_term_ups(up)) {
          auto [up_flip, coeff] = term_action(up);
          if ((up_flip & dn) == 0) { // tJ constraint
            int64_t idx_out = basis_out.index_r(up_flip, dn);
            vec_out[idx_out] += coeff * vec_in[idx_in];
          } // tJ constraint
        } // non-zero term dns
      } // if ((upspins & flipmask) == 0)
    } // non-zero-term ups

    ++idx_dn;
  } // for(const bit_t& upspins : my_upspins_)
}

} // namespace xdiag::basis::tj_distributed
