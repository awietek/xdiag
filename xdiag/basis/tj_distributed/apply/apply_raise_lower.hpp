#pragma once

#include <string>

#include <xdiag/basis/tj_distributed/apply/generic_term_dns.hpp>
#include <xdiag/basis/tj_distributed/apply/generic_term_ups.hpp>
#include <xdiag/bits/bitops.hpp>

namespace xdiag::basis::tj_distributed {

template <typename bit_t, typename coeff_t, class BasisIn, class BasisOut,
          class Fill>
void apply_raise_lower(Coupling const &cpl, Op const &op, BasisIn &&basis_in,
                       BasisOut &&basis_out, Fill &&fill) {
  coeff_t c = cpl.scalar().as<coeff_t>();
  std::string type = op.type();
  int64_t s = op[0];

  bit_t site_mask = (bit_t)1 << s;
  bit_t fermi_mask = site_mask - 1;

  // Raising operators
  if ((type == "Cdagup") || (type == "Cdagdn")) {

    auto term_action = [&](bit_t spins) -> std::pair<bit_t, coeff_t> {
      bool fermi = bits::popcnt(spins & fermi_mask) & 1;
      return {spins ^ site_mask, fermi ? -c : c};
    };

    if (type == "Cdagup") {
      auto non_zero_term = [&](bit_t const &ups) -> bool {
        return (ups & site_mask) == 0;
      };
      tj_distributed::generic_term_ups<bit_t, coeff_t>(
          basis_in, basis_out, non_zero_term, term_action, fill);
    } else if (type == "Cdagdn") {
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
  } else if ((type == "Cup") || (type == "Cdn")) {

    auto term_action = [&](bit_t spins) -> std::pair<bit_t, coeff_t> {
      bool fermi = bits::popcnt(spins & fermi_mask) & 1;
      return {spins ^ site_mask, fermi ? -c : c};
    };

    if (type == "Cup") {
      auto non_zero_term = [&](bit_t const &spins) -> bool {
        return (spins & site_mask);
      };
      tj_distributed::generic_term_ups<bit_t, coeff_t>(
          basis_in, basis_out, non_zero_term, term_action, fill);
    } else if (type == "Cdn") {
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
} // namespace xdiag::basis::tj_distributed
