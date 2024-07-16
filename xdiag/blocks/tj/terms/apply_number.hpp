#pragma once

#include <xdiag/blocks/tj/terms/generic_term_diag.hpp>
#include <xdiag/common.hpp>
#include <xdiag/operators/op.hpp>

namespace xdiag::tj {

template <typename bit_t, typename coeff_t, bool symmetric, class Basis,
          class Fill>
void apply_number(Op const &op, Basis &&basis, Fill &&fill) {
  assert(op.coupling_defined());
  assert(op.type_defined());
  assert(op.size() == 1);

  std::string type = op.type();
  assert((type == "NUMBERUP") || (type == "NUMBERDN"));

  Coupling cpl = op.coupling();
  assert(cpl.isexplicit() && !cpl.ismatrix());
  coeff_t mu = cpl.as<coeff_t>();

  int64_t s = op[0];
  bit_t mask = (bit_t)1 << s;

  if (type == "NUMBERUP") {
    auto term_action = [&](bit_t up, bit_t dn) {
      (void)dn;
      return (up & mask) ? mu : 0.;
    };
    tj::generic_term_diag<bit_t, coeff_t, symmetric>(basis, term_action, fill);
  } else if (type == "NUMBERDN") {
    auto term_action = [&](bit_t up, bit_t dn) {
      (void)up;
      return (dn & mask) ? mu : 0.;
    };
    tj::generic_term_diag<bit_t, coeff_t, symmetric>(basis, term_action, fill);
  }
}

} // namespace xdiag::tj
