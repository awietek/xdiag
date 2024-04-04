#pragma once

#include <vector>

#include <xdiag/bits/bitops.hpp>

namespace xdiag::tj_distributed {

template <typename bit_t, typename coeff_t, class BasisIn, class BasisOut,
          class NonZeroTermUps, class NonZeroTermDns, class TermAction>
void generic_term_ups(BasisIn &&basis_in, BasisOut &&basis_out,
                      NonZeroTermUps &&non_zero_term_ups,
                      NonZeroTermDns &&non_zero_term_dns,
                      TermAction &&term_action, const coeff_t *vec_in,
                      coeff_t *vec_out) try {

  int64_t n_sites = basis_in.n_sites();
  assert(n_sites == basis_out.n_sites());
  bit_t sitesmask = ((bit_t)1 << n_sites) - 1;

  int64_t n_up_in = basis_in.n_up();
  int64_t n_dn_in = basis_in.n_dn();
  int64_t n_up_configurations_in =
      combinatorics::binomial(n_sites - n_dn_in, n_up_in);

  int64_t n_up_out = basis_out.n_up();
  int64_t n_dn_out = basis_out.n_dn();
  int64_t n_up_configurations_out =
      combinatorics::binomial(n_sites - n_dn_out, n_up_out);

  // Loop over all configurations
  int64_t idx_dn = 0;
  for (bit_t dn : basis_in.my_dns()) {

    if (non_zero_term_dns(dn)) {
      bit_t not_dn = (~dn) & sitesmask;
      int64_t dn_offset_in = idx_dn * n_up_configurations_in;
      int64_t dn_offset_out = idx_dn * n_up_configurations_out;

      for (int64_t idx_in = dn_offset_in;
           idx_in < dn_offset_in + n_up_configurations_in; ++idx_in) {

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
} catch (...) {
  XDiagRethrow("Unable to apply generic term on up spins");
}

} // namespace xdiag::tj_distributed
