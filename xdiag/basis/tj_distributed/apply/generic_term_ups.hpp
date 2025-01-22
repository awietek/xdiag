#pragma once

#include <vector>

#include <xdiag/bits/bitops.hpp>

namespace xdiag::basis::tj_distributed {

template <typename bit_t, typename coeff_t, class BasisIn, class BasisOut,
          class NonZeroTermUps, class NonZeroTermDns, class TermAction>
void generic_term_ups(BasisIn &&basis_in, BasisOut &&basis_out,
                      NonZeroTermUps &&non_zero_term_ups,
                      NonZeroTermDns &&non_zero_term_dns,
                      TermAction &&term_action, const coeff_t *vec_in,
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
      int64_t dn_offset_out = idx_dn * nup_configurations_out;

      for (int64_t idx_in = dn_offset_in;
           idx_in < dn_offset_in + nup_configurations_in; ++idx_in) {

        bit_t up = basis_in.my_ups_for_dns_storage(idx_in);

        // Check if hopping is possible
        if (non_zero_term_ups(up)) {
          auto [up_flip, coeff] = term_action(up);
          if ((up_flip & dn) == 0) { // tJ constraint
            bit_t upc_flip = bits::extract(up_flip, not_dn);
            int64_t idx_upc_flip = basis_out.index_upcs(upc_flip);
            int64_t idx_out = dn_offset_out + idx_upc_flip;

            vec_out[idx_out] += coeff * vec_in[idx_in];
          } // tJ constraint
        }   // non-zero term dns
      }     // if ((upspins & flipmask) == 0)
    }       // non-zero-term ups

    ++idx_dn;
  } // for(const bit_t& upspins : my_upspins_)
}

} // namespace xdiag::basis::tj_distributed
