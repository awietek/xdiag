#pragma once

#include <xdiag/bits/bitops.hpp>
#include <xdiag/common.hpp>
#include <xdiag/operators/op.hpp>

#include <xdiag/basis/electron/apply/generic_term_dns.hpp>
#include <xdiag/basis/electron/apply/generic_term_ups.hpp>

namespace xdiag::basis::electron {

template <typename bit_t, typename coeff_t, bool symmetric, class Basis,
          class Fill>
void apply_hopping(Coupling const &cpl, Op const &op, Basis &&basis,
                   Fill &&fill) try {

  coeff_t t = cpl.scalar().as<coeff_t>();
  int64_t s1 = op[0];
  int64_t s2 = op[1];
  bit_t flipmask = ((bit_t)1 << s1) | ((bit_t)1 << s2);
  int64_t l = std::min(s1, s2);
  int64_t u = std::max(s1, s2);
  bit_t fermimask = (((bit_t)1 << (u - l - 1)) - 1) << (l + 1);

  auto non_zero_term = [&flipmask](bit_t const &spins) -> bool {
    return bits::popcnt(spins & flipmask) & 1;
  };

  auto term_action = [&](bit_t spins) -> std::pair<bit_t, coeff_t> {
    bool fermi = bits::popcnt(spins & fermimask) & 1;
    spins ^= flipmask;
    if constexpr (isreal<coeff_t>()) {
      return {spins, fermi ? t : -t};
    } else {
      coeff_t tt = (bits::gbit(spins, s1)) ? t : xdiag::conj(t);
      return {spins, fermi ? tt : -tt};
    }
  };

  std::string type = op.type();
  if (type == "Hopup") {
    electron::generic_term_ups<bit_t, coeff_t, symmetric>(
        basis, basis, non_zero_term, term_action, fill);
  } else if (type == "Hopdn") {
    electron::generic_term_dns<bit_t, coeff_t, symmetric, false>(
        basis, basis, non_zero_term, term_action, fill);
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

} // namespace xdiag::basis::electron
