#pragma once

#include <xdiag/bits/bitops.h>
#include <xdiag/common.h>
#include <xdiag/operators/bond.h>

#include <xdiag/blocks/electron/terms/generic_term_ups.h>
#include <xdiag/blocks/electron/terms/generic_term_dns.h>

namespace xdiag::electron {

template <typename bit_t, typename coeff_t, bool symmetric, class Basis,
          class Fill>
void apply_hopping(Bond const &bond, Basis &&basis, Fill &&fill) {
  assert(bond.coupling_defined());
  assert(bond.type_defined());
  assert(bond.size() == 2);
  assert(bond.sites_disjoint());

  std::string type = bond.type();
  assert((type == "HOPUP") || (type == "HOPDN"));

  int64_t s1 = bond[0];
  int64_t s2 = bond[1];
  bit_t flipmask = ((bit_t)1 << s1) | ((bit_t)1 << s2);
  int64_t l = std::min(s1, s2);
  int64_t u = std::max(s1, s2);
  bit_t fermimask = (((bit_t)1 << (u - l - 1)) - 1) << (l + 1);
  coeff_t t = bond.coupling<coeff_t>();

  auto non_zero_term = [&flipmask](bit_t const &spins) -> bool {
    return bits::popcnt(spins & flipmask) & 1;
  };

  auto term_action = [&](bit_t spins) -> std::pair<bit_t, coeff_t> {
    bool fermi = bits::popcnt(spins & fermimask) & 1;
    spins ^= flipmask;
    if constexpr (iscomplex<coeff_t>()) {
      coeff_t tt = (bits::gbit(spins, s1)) ? t : xdiag::conj(t);
      return {spins, fermi ? tt : -tt};
    } else {
      return {spins, fermi ? t : -t};
    }
  };

  if (type == "HOPUP") {
    electron::generic_term_ups<bit_t, coeff_t, symmetric>(
        basis, basis, non_zero_term, term_action, fill);
  } else if (type == "HOPDN") {
    electron::generic_term_dns<bit_t, coeff_t, symmetric, false>(
        basis, basis, non_zero_term, term_action, fill);
  }
  
}

} // namespace xdiag::electron
