#pragma once

#include <xdiag/basis/electron_distributed/apply/generic_term_diag.hpp>
#include <xdiag/bits/bitops.hpp>
#include <xdiag/common.hpp>
#include <xdiag/operators/op.hpp>

namespace xdiag::basis::electron_distributed {

template <typename bit_t, typename coeff_t, class Basis>
void apply_szsz(Coupling const &cpl, Op const &op, Basis &&basis,
                const coeff_t *vec_in, coeff_t *vec_out) {
  coeff_t J = cpl.scalar().as<coeff_t>();
  int64_t s1 = op[0];
  int64_t s2 = op[1];
  bit_t mask = ((bit_t)1 << s1) | ((bit_t)1 << s2);

  coeff_t val_same = J / 4.;
  coeff_t val_diff = -J / 4.;

  auto term_action = [&](bit_t up, bit_t dn) -> coeff_t {
    // one site either empty or doubly occupied
    if (((up ^ dn) & mask) != mask) {
      return 0.;
    } else {
      return bits::popcnt(up & mask) == 1 ? val_diff : val_same;
    }
  };

  electron_distributed::generic_term_diag<bit_t, coeff_t>(basis, term_action,
                                                          vec_in, vec_out);
}

} // namespace xdiag::basis::electron_distributed
