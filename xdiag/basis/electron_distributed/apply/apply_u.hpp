#pragma once

#include <xdiag/basis/electron_distributed/apply/generic_term_diag.hpp>
#include <xdiag/common.hpp>
#include <xdiag/operators/op.hpp>

namespace xdiag::basis::electron_distributed {

template <typename bit_t, typename coeff_t, class basis_t>
void apply_u(Coupling const &cpl, basis_t const &basis, const coeff_t *vec_in,
             coeff_t *vec_out) {
  coeff_t U = cpl.scalar().as<coeff_t>();
  auto term_action = [&](bit_t up, bit_t dn) {
    return U * (double)bits::popcnt(up & dn);
  };
  electron_distributed::generic_term_diag<bit_t, coeff_t>(basis, term_action,
                                                          vec_in, vec_out);
}

} // namespace xdiag::basis::electron_distributed
