#pragma once

#include <xdiag/basis/tj/apply/generic_term_mixed.hpp>
#include <xdiag/bits/bitops.hpp>
#include <xdiag/common.hpp>
#include <xdiag/operators/op.hpp>

namespace xdiag::basis::tj {

template <typename coeff_t, bool symmetric, class basis_t, class fill_f>
void apply_exchange(Coupling const &cpl, Op const &op, basis_t const &basis,
                    fill_f fill) try {
  using bit_t = typename basis_t::bit_t;

  coeff_t J = cpl.scalar().as<coeff_t>();
  coeff_t Jhalf = J / 2.;
  coeff_t Jhalf_conj = conj(Jhalf);

  int64_t s1 = op[0];
  int64_t s2 = op[1];

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
  auto term_actionups = [&](bit_t up) -> std::pair<bit_t, coeff_t> {
    bool fermi_up = bits::popcnt(up & fermimask) & 1;
    bit_t up_flip = up ^ flipmask;
    if constexpr (isreal<coeff_t>()) {
      return {up_flip, fermi_up ? Jhalf : -Jhalf};
    } else {
      if (bits::gbit(up, s1)) {
        return {up_flip, fermi_up ? Jhalf : -Jhalf};
      } else {
        return {up_flip, fermi_up ? Jhalf_conj : -Jhalf_conj};
      }
    }
  };
  auto term_actiondns = [&](bit_t dn) -> std::pair<bit_t, coeff_t> {
    bool fermi_dn = bits::popcnt(dn & fermimask) & 1;
    return {dn ^ flipmask, fermi_dn ? -1.0 : 1.0};
  };
  generic_term_mixed<bit_t, coeff_t, symmetric>(
      basis, basis, non_zero_term_ups, non_zero_term_dns, term_actionups,
      term_actiondns, fill);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

} // namespace xdiag::basis::tj
