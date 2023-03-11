#pragma once

#include <hydra/bitops/bitops.h>
#include <hydra/common.h>
#include <hydra/operators/bond.h>

#include <hydra/blocks/tj/terms/generic_term_diag.h>

namespace hydra::electron {

template <typename bit_t, typename coeff_t, bool symmetric, class Indexing,
          class Fill>
void apply_number(Bond const &bond, Indexing &&indexing, Fill &&fill) {
  assert(bond.coupling_defined());
  assert(bond.type_defined());
  assert(bond.size() == 1);

  std::string type = bond.type();
  assert((type == "NUMBERUP") || (type == "NUMBERDN"));

  coeff_t mu = bond.coupling<coeff_t>();
  int s = bond.site(0);
  bit_t mask = (bit_t)1 << s;

  if (type == "NUMBERUP") {
    auto term_action = [&](bit_t up, bit_t dn) {
      return (up & mask) ? mu : 0.;
    };
    tj::generic_term_diag(indexing, term_action, fill);
  } else if (type == "NUMBERDN") {
    auto term_action = [&](bit_t up, bit_t dn) {
      return (dn & mask) ? mu : 0.;
    };
    tj::generic_term_diag(indexing, term_action, fill);
  }
}

} // namespace hydra::electron
