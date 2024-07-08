#pragma once

#include <xdiag/blocks/tj/terms/generic_term_diag.hpp>
#include <xdiag/common.hpp>
#include <xdiag/operators/bond.hpp>

namespace xdiag::tj_distributed {

template <typename bit_t, typename coeff_t, class Basis>
void apply_number(Bond const &bond, Basis &&basis, const coeff_t *vec_in,
                  coeff_t *vec_out) {
  assert(bond.coupling_defined());
  assert(bond.type_defined());
  assert(bond.size() == 1);

  std::string type = bond.type();
  assert((type == "NUMBERUP") || (type == "NUMBERDN"));

  coeff_t mu = bond.coupling<coeff_t>();
  int64_t s = bond.site(0);
  bit_t mask = (bit_t)1 << s;

  if (type == "NUMBERUP") {
    auto term_action = [&](bit_t up, bit_t dn) {
      (void)dn;
      return (up & mask) ? mu : 0.;
    };
    tj_distributed::generic_term_diag<bit_t, coeff_t>(basis, term_action,
                                                      vec_in, vec_out);
  } else if (type == "NUMBERDN") {
    auto term_action = [&](bit_t up, bit_t dn) {
      (void)up;
      return (dn & mask) ? mu : 0.;
    };
    tj_distributed::generic_term_diag<bit_t, coeff_t>(basis, term_action,
                                                      vec_in, vec_out);
  }
}

} // namespace xdiag::tj_distributed
