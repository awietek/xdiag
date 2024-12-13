#pragma once

#include <xdiag/basis/tj/apply/generic_term_dns.hpp>
#include <xdiag/basis/tj/apply/generic_term_ups.hpp>
#include <xdiag/bits/bitops.hpp>
#include <xdiag/common.hpp>
#include <xdiag/operators/op.hpp>

namespace xdiag::basis::tj {

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

  auto term_action = [&](bit_t spins) -> std::pair<bit_t, coeff_t> {
    bool fermi = bits::popcnt(spins & fermimask) & 1;
    spins ^= flipmask;
    if constexpr (iscomplex<coeff_t>()) {
      coeff_t tt = (bits::gbit(spins, s1)) ? t : conj(t);
      return {spins, fermi ? tt : -tt};
    } else {
      return {spins, fermi ? t : -t};
    }
  };

  std::string type = op.type();
  if (type == "Hopup") {
    auto non_zero_term = [&flipmask](bit_t const &ups) -> bool {
      return bits::popcnt(ups & flipmask) & 1;
    };
    tj::generic_term_ups<bit_t, coeff_t, symmetric>(basis, basis, non_zero_term,
                                                    term_action, fill);
  } else if (type == "Hopdn") {
    auto non_zero_term_ups = [&flipmask](bit_t const &ups) -> bool {
      return (ups & flipmask) == 0;
    };
    auto non_zero_term_dns = [&flipmask](bit_t const &dns) -> bool {
      return bits::popcnt(dns & flipmask) & 1;
    };
    tj::generic_term_dns<bit_t, coeff_t, symmetric, false>(
        basis, basis, non_zero_term_ups, non_zero_term_dns, term_action, fill);
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

} // namespace xdiag::basis::tj
