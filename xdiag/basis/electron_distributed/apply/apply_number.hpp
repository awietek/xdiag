#pragma once

#include <xdiag/basis/electron_distributed/apply/generic_term_diag.hpp>
#include <xdiag/common.hpp>
#include <xdiag/operators/op.hpp>

namespace xdiag::basis::electron_distributed {

template <typename bit_t, typename coeff_t, class Basis>
void apply_number(Coupling const &cpl, Op const &op, Basis &&basis,
                  const coeff_t *vec_in, coeff_t *vec_out) {
  coeff_t mu = cpl.scalar().as<coeff_t>();
  int64_t s = op[0];
  std::string type = op.type();

  bit_t mask = (bit_t)1 << s;
  if (type == "Nup") {
    auto term_action = [&](bit_t up, bit_t dn) {
      (void)dn;
      return (up & mask) ? mu : 0.;
    };
    electron_distributed::generic_term_diag<bit_t, coeff_t>(basis, term_action,
                                                            vec_in, vec_out);
  } else if (type == "Ndn") {
    auto term_action = [&](bit_t up, bit_t dn) {
      (void)up;
      return (dn & mask) ? mu : 0.;
    };
    electron_distributed::generic_term_diag<bit_t, coeff_t>(basis, term_action,
                                                            vec_in, vec_out);
  }
}

} // namespace xdiag::basis::electron_distributed
