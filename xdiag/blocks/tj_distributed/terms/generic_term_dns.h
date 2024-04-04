#pragma once

#include <functional>
#include <xdiag/bits/bitops.h>

namespace xdiag::tj_distributed {

template <typename bit_t, typename coeff_t, bool fermi_ups, class BasisIn,
          class BasisOut, class NonZeroTermUps, class NonZeroTermDns,
          class TermAction>
void generic_term_dns(BasisIn &&basis_in, BasisOut &&basis_out,
                      NonZeroTermUps &&non_zero_term_ups,
                      NonZeroTermDns &&non_zero_term_dns,
                      TermAction &&term_action, const coeff_t *vec_in,
                      coeff_t *vec_out) try {

  int64_t n_sites = basis_in.n_sites();
  assert(n_sites == basis_out.n_sites());
  bit_t sitesmask = ((bit_t)1 << n_sites) - 1;

  int64_t n_up_in = basis_in.n_up();
  int64_t n_dn_in = basis_in.n_dn();
  int64_t n_dn_configurations_in =
      combinatorics::binomial(n_sites - n_up_in, n_dn_in);

  int64_t n_up_out = basis_out.n_up();
  int64_t n_dn_out = basis_out.n_dn();
  int64_t n_dn_configurations_out =
      combinatorics::binomial(n_sites - n_up_out, n_dn_out);

  // Loop over all configurations
  int64_t idx_up = 0;
  for (bit_t up : basis_in.my_ups()) {

    if (non_zero_term_ups(up)) {
      bit_t not_up = (~up) & sitesmask;
      int64_t up_offset_in = idx_up * n_dn_configurations_in;
      int64_t up_offset_out = idx_up * n_dn_configurations_out;

      for (int64_t idx_in = up_offset_in;
           idx_in < up_offset_in + n_dn_configurations_in; ++idx_in) {

        bit_t dn = basis_in.my_dns_for_ups_storage(idx_in);

        // Check if hopping is possible
        if (non_zero_term_dns(dn)) {
          auto [dn_flip, coeff] = term_action(dn);
          if ((dn_flip & up) == 0) { // tJ constraint
            bit_t dnc_flip = bits::extract(dn_flip, not_up);
            int64_t idx_dnc_flip = basis_out.index_dncs(dnc_flip);
            int64_t idx_out = up_offset_out + idx_dnc_flip;

            if constexpr (fermi_ups) { // if odd number of ups -> factor -1
              bool fermi_up = (bool)(bits::popcnt(up) & 1);
              vec_out[idx_out] += (fermi_up ? -coeff : coeff) * vec_in[idx_in];
            } else {
              vec_out[idx_out] += coeff * vec_in[idx_in];
            }
          } // tJ constraint
        }   // non-zero term dns
      }     // if ((upspins & flipmask) == 0)
    }       // non-zero-term ups

    ++idx_up;
  } // for(const bit_t& upspins : my_upspins_)
} catch (...) {
  XDiagRethrow("Unable to apply generic term on dn spins");
}

} // namespace xdiag::tj_distributed
