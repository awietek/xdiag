// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <functional>

#include <xdiag/bits/bitops.hpp>

namespace xdiag::basis::electron_distributed {

template <typename coeff_t, bool fermi_ups, class basis_t,
          class non_zero_term_ups_f, class non_zero_term_dns_f,
          class term_action_f>
void generic_term_dns(basis_t const &basis_in, basis_t const &basis_out,
                      non_zero_term_ups_f non_zero_term_ups,
                      non_zero_term_dns_f non_zero_term_dns,
                      term_action_f term_action, const coeff_t *vec_in,
                      coeff_t *vec_out) {
  using bit_t = typename basis_t::bit_t;

  int64_t nsites = basis_in.nsites();
  assert(nsites == basis_out.nsites());
  bit_t sitesmask = ((bit_t)1 << nsites) - 1;

  int64_t nup_in = basis_in.nup();
  int64_t ndn_in = basis_in.ndn();
  int64_t ndn_configurations_in = combinatorics::binomial(nsites, ndn_in);

  int64_t nup_out = basis_out.nup();
  int64_t ndn_out = basis_out.ndn();
  int64_t ndn_configurations_out = combinatorics::binomial(nsites, ndn_out);

  // Loop over all configurations
  int64_t idx_up = 0;
  for (bit_t up : basis_in.my_ups()) {

    if (non_zero_term_ups(up)) {
      int64_t up_offset_in = idx_up * ndn_configurations_in;
      int64_t up_offset_out = idx_up * ndn_configurations_out;

      int64_t idx_in = up_offset_in;
      for (bit_t dn : basis_in.all_dns()) {

        // Check if hopping is possible
        if (non_zero_term_dns(dn)) {
          auto [dn_flip, coeff] = term_action(dn);
          int64_t idx_dn_flip = basis_out.index_dns(dn_flip);
          int64_t idx_out = up_offset_out + idx_dn_flip;

          if constexpr (fermi_ups) { // if odd number of ups -> factor -1
            bool fermi_up = (bool)(bits::popcnt(up) & 1);
            vec_out[idx_out] += (fermi_up ? -coeff : coeff) * vec_in[idx_in];
          } else {
            vec_out[idx_out] += coeff * vec_in[idx_in];
          }
        } // non-zero term dns
        ++idx_in;
      }
    } // non-zero-term ups
    ++idx_up;
  }
}

} // namespace xdiag::basis::electron_distributed
