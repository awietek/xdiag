#pragma once

#include <hydra/bitops/bitops.h>
#include <hydra/common.h>
#include <hydra/operators/bond.h>

#include <hydra/blocks/tj/terms/generic_term_mixed.h>

namespace hydra::tj {

template <typename bit_t, typename coeff_t, bool symmetric, class Indexing,
          class Filler>
void apply_exchange(Bond const &bond, Indexing &&indexing, Filler &&fill) {
  assert(bond.coupling_defined());
  assert(bond.type_defined());
  assert(bond.size() == 2);
  assert(bond.sites_disjoint());
  std::string type = bond.type();
  assert(type == "EXCHANGE");

  int s1 = bond[0];
  int s2 = bond[1];
  coeff_t J = bond.coupling<coeff_t>();
  coeff_t Jhalf = J / 2.;
  coeff_t Jhalf_conj = hydra::conj(J) / 2.;

  // Prepare bitmasks
  bit_t flipmask = ((bit_t)1 << s1) | ((bit_t)1 << s2);
  // bit_t sitesmask = ((bit_t)1 << n_sites) - 1;
  // int l = std::min(s1, s2);
  // int u = std::max(s1, s2);
  // bit_t fermimask = (((bit_t)1 << (u - l - 1)) - 1) << (l + 1);

  auto non_zero_term_ups = [&](bit_t up) {
    return bitops::popcnt(up & flipmask) == 1;
  };
  auto non_zero_term_dns = [&](bit_t dn) {
    return bitops::popcnt(dn & flipmask) == 1;
  };
  auto term_action_ups = [&](bit_t up) { return up ^ flipmask; };
  auto term_action_dns = [&](bit_t dn) { return dn ^ flipmask; };
  tj::generic_term_mixed(indexing, indexing, non_zero_term_ups,
                         non_zero_term_dns, term_action_ups, term_action_dns,
                         fill);

} // namespace hydra::electron
