#pragma once

#include <xdiag/basis/tj/apply/generic_term_diag.hpp>
#include <xdiag/common.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/bits/gbit.hpp>

namespace xdiag::basis::tj {

template <typename bit_t, typename coeff_t, bool symmetric, class Basis,
          class Fill>
void apply_number_number(Coupling const &cpl, Op const &op, Basis &&basis,
                         Fill &&fill) try {
  using bits::gbit;
  coeff_t mu = cpl.scalar().as<coeff_t>();
  int64_t s1 = op[0];
  int64_t s2 = op[1];

  auto term_action = [&](bit_t up, bit_t dn) {
    int n1 = gbit(up, s1) + gbit(dn, s1);
    int n2 = gbit(up, s2) + gbit(dn, s2);
    return mu * (coeff_t)(n1 * n2);
  };
  tj::generic_term_diag<bit_t, coeff_t, symmetric>(basis, term_action, fill);

} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

} // namespace xdiag::basis::tj
