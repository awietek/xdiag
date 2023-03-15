#pragma once

#include <string>

#include <hydra/bitops/bitops.h>
#include <hydra/blocks/tj/terms/generic_term_dns.h>
#include <hydra/blocks/tj/terms/generic_term_ups.h>

namespace hydra::tj {

template <typename bit_t, typename coeff_t, bool symmetric, class IndexingIn,
          class IndexingOut, class Fill>
void apply_raise_lower(Bond const &bond, IndexingIn &&indexing_in,
                       IndexingOut &&indexing_out, Fill &&fill) {
  assert(bond.coupling_defined());
  assert(bond.type_defined());
  assert(bond.size() == 1);

  std::string type = bond.type();
  assert((type == "CDAGUP") || (type == "CDAGDN") || (type == "CUP") ||
         (type == "CDN"));
  
  int s = bond[0];
  bit_t site_mask = (bit_t)1 << s;
  bit_t fermi_mask = site_mask - 1;
  coeff_t c = bond.coupling<coeff_t>();

  // Raising operators
  if ((type == "CDAGUP") || (type == "CDAGDN")) {

    auto term_action = [&](bit_t spins) -> std::pair<bit_t, coeff_t> {
      bool fermi = bitops::popcnt(spins & fermi_mask) & 1;
      return {spins ^ site_mask, fermi ? -c : c};
    };

    if (type == "CDAGUP") {
      auto non_zero_term = [&](bit_t const &ups) -> bool {
        return (ups & site_mask) == 0;
      };
      tj::generic_term_ups<bit_t, coeff_t, symmetric>(
          indexing_in, indexing_out, non_zero_term, term_action, fill);
    } else if (type == "CDAGDN") {
      auto non_zero_term_ups = [&](bit_t const &ups) -> bool {
        return (ups & site_mask) == 0;
      };
      auto non_zero_term_dns = [&](bit_t const &dns) -> bool {
        return (dns & site_mask) == 0;
      };
      tj::generic_term_dns<bit_t, coeff_t, symmetric, true>(
          indexing_in, indexing_out, non_zero_term_ups, non_zero_term_dns,
          term_action, fill);
    }

    // Lowering operators
  } else if ((type == "CUP") || (type == "CDN")) {

    auto term_action = [&](bit_t spins) -> std::pair<bit_t, coeff_t> {
      bool fermi = bitops::popcnt(spins & fermi_mask) & 1;
      return {spins ^ site_mask, fermi ? -c : c};
    };

    if (type == "CUP") {
      auto non_zero_term = [&](bit_t const &spins) -> bool {
        return (spins & site_mask);
      };
      tj::generic_term_ups<bit_t, coeff_t, symmetric>(
          indexing_in, indexing_out, non_zero_term, term_action, fill);
    } else if (type == "CDN") {
      auto non_zero_term_ups = [&](bit_t const &ups) -> bool {
        return (ups & site_mask) == 0;
      };
      auto non_zero_term_dns = [&](bit_t const &dns) -> bool {
        return (dns & site_mask);
      };
      tj::generic_term_dns<bit_t, coeff_t, symmetric, true>(
          indexing_in, indexing_out, non_zero_term_ups, non_zero_term_dns,
          term_action, fill);
    }
  }
}
} // namespace hydra::tj
