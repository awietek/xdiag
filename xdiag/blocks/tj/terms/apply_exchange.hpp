#pragma once

#include <xdiag/bits/bitops.hpp>
#include <xdiag/blocks/tj/terms/generic_term_mixed.hpp>
#include <xdiag/common.hpp>
#include <xdiag/operators/bond.hpp>

namespace xdiag::tj {

template <typename bit_t, typename coeff_t, bool symmetric, class Basis,
          class Filler>
void apply_exchange(Bond const &bond, Basis &&basis, Filler &&fill) {
  assert(bond.coupling_defined());
  assert(bond.type_defined());
  assert(bond.size() == 2);
  assert(bond.sites_disjoint());
  std::string type = bond.type();
  assert(type == "EXCHANGE");

  int64_t s1 = bond[0];
  int64_t s2 = bond[1];
  coeff_t J = bond.coupling<coeff_t>();
  coeff_t Jhalf = J / 2.;
  coeff_t Jhalf_conj = conj(Jhalf);

  // Prepare bitmasks
  bit_t flipmask = ((bit_t)1 << s1) | ((bit_t)1 << s2);
  int64_t l = std::min(s1, s2);
  int64_t u = std::max(s1, s2);
  bit_t fermimask = (((bit_t)1 << (u - l - 1)) - 1) << (l + 1);

  auto non_zero_term_ups = [&](bit_t up) {
    return bits::popcnt(up & flipmask) == 1;
  };
  auto non_zero_term_dns = [&](bit_t dn) {
    return bits::popcnt(dn & flipmask) == 1;
  };
  auto term_action_ups = [&](bit_t up) -> std::pair<bit_t, coeff_t> {
    bool fermi_up = bits::popcnt(up & fermimask) & 1;
    bit_t up_flip = up ^ flipmask;
    if constexpr (iscomplex<coeff_t>()) {
      if (bits::gbit(up, s1)) {
        // Log("a {}", Jhalf);
        return {up_flip, fermi_up ? Jhalf : -Jhalf};
      } else {
        // Log("b {}", Jhalf_conj);
        return {up_flip, fermi_up ? Jhalf_conj : -Jhalf_conj};
      }
    } else {
      return {up_flip, fermi_up ? Jhalf : -Jhalf};
    }
  };
  auto term_action_dns = [&](bit_t dn) -> std::pair<bit_t, coeff_t> {
    bool fermi_dn = bits::popcnt(dn & fermimask) & 1;
    return {dn ^ flipmask, fermi_dn ? -1.0 : 1.0};
  };
  tj::generic_term_mixed<bit_t, coeff_t, symmetric>(
      basis, basis, non_zero_term_ups, non_zero_term_dns, term_action_ups,
      term_action_dns, fill);
}

} // namespace xdiag::tj
