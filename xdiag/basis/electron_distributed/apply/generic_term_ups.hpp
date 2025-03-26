#pragma once

#include <vector>

#include <xdiag/bits/bitops.hpp>

namespace xdiag::basis::electron_distributed {

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
  int64_t nup_configurations_in = combinatorics::binomial(nsites, nup_in);

  int64_t nup_out = basis_out.nup();
  int64_t ndn_out = basis_out.ndn();
  int64_t nup_configurations_out = combinatorics::binomial(nsites, nup_out);

  // Loop over all configurations
  int64_t idx_dn = 0;
  for (bit_t dn : basis_in.my_dns()) {

    if (non_zero_term_dns(dn)) {
      int64_t dn_offset_in = idx_dn * nup_configurations_in;
      int64_t dn_offset_out = idx_dn * nup_configurations_out;

      int64_t idx_in = dn_offset_in;
      for (bit_t up : basis_in.all_ups()) {

        // Check if hopping is possible
        if (non_zero_term_ups(up)) {
          auto [up_flip, coeff] = term_action(up);
          int64_t idx_up_flip = basis_out.index_ups(up_flip);
          int64_t idx_out = dn_offset_out + idx_up_flip;
          vec_out[idx_out] += coeff * vec_in[idx_in];
        } // non-zero term dns

        ++idx_in;
      } // for (bit_t up : basis_in.all_ups())
    } // non-zero-term ups

    ++idx_dn;
  } // for (bit_t dn : basis_in.my_dns())
}

} // namespace xdiag::basis::electron_distributed
