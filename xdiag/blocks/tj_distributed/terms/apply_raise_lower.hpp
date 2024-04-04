#pragma once

#include <string>

#include <xdiag/bits/bitops.hpp>
#include <xdiag/blocks/tj_distributed/terms/generic_term_dns.hpp>
#include <xdiag/blocks/tj_distributed/terms/generic_term_ups.hpp>

namespace xdiag::tj_distributed {

template <typename bit_t, typename coeff_t, class BasisIn, class BasisOut,
          class Fill>
void apply_raise_lower(Bond const &bond, BasisIn &&basis_in,
                       BasisOut &&basis_out, Fill &&fill) {
  assert(bond.coupling_defined());
  assert(bond.type_defined());
  assert(bond.size() == 1);

  std::string type = bond.type();
  assert((type == "CDAGUP") || (type == "CDAGDN") || (type == "CUP") ||
         (type == "CDN"));

  int64_t s = bond[0];
  bit_t site_mask = (bit_t)1 << s;
  bit_t fermi_mask = site_mask - 1;
  coeff_t c = bond.coupling<coeff_t>();

  // Raising operators
  if ((type == "CDAGUP") || (type == "CDAGDN")) {

    auto term_action = [&](bit_t spins) -> std::pair<bit_t, coeff_t> {
      bool fermi = bits::popcnt(spins & fermi_mask) & 1;
      return {spins ^ site_mask, fermi ? -c : c};
    };

    if (type == "CDAGUP") {
      auto non_zero_term = [&](bit_t const &ups) -> bool {
        return (ups & site_mask) == 0;
      };
      tj_distributed::generic_term_ups<bit_t, coeff_t>(
          basis_in, basis_out, non_zero_term, term_action, fill);
    } else if (type == "CDAGDN") {
      auto non_zero_term_ups = [&](bit_t const &ups) -> bool {
        return (ups & site_mask) == 0;
      };
      auto non_zero_term_dns = [&](bit_t const &dns) -> bool {
        return (dns & site_mask) == 0;
      };
      tj_distributed::generic_term_dns<bit_t, coeff_t, true>(
          basis_in, basis_out, non_zero_term_ups, non_zero_term_dns,
          term_action, fill);
    }

    // Lowering operators
  } else if ((type == "CUP") || (type == "CDN")) {

    auto term_action = [&](bit_t spins) -> std::pair<bit_t, coeff_t> {
      bool fermi = bits::popcnt(spins & fermi_mask) & 1;
      return {spins ^ site_mask, fermi ? -c : c};
    };

    if (type == "CUP") {
      auto non_zero_term = [&](bit_t const &spins) -> bool {
        return (spins & site_mask);
      };
      tj_distributed::generic_term_ups<bit_t, coeff_t>(
          basis_in, basis_out, non_zero_term, term_action, fill);
    } else if (type == "CDN") {
      auto non_zero_term_ups = [&](bit_t const &ups) -> bool {
        return (ups & site_mask) == 0;
      };
      auto non_zero_term_dns = [&](bit_t const &dns) -> bool {
        return (dns & site_mask);
      };
      tj_distributed::generic_term_dns<bit_t, coeff_t, true>(
          basis_in, basis_out, non_zero_term_ups, non_zero_term_dns,
          term_action, fill);
    }
  }
}
} // namespace xdiag::tj_distributed
