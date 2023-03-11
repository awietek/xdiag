#pragma once

#include <hydra/bitops/bitops.h>
#include <hydra/common.h>
#include <hydra/operators/bond.h>

#include <hydra/blocks/tj/terms/generic_term_diag.h>

namespace hydra::electron {

template <typename bit_t, typename coeff_t, bool symmetric, class Indexing,
          class Fill>
void apply_ising(Bond const &bond, Indexing &&indexing, Fill &&fill) {
  assert(bond.coupling_defined());
  assert(bond.type_defined());
  assert(bond.size() == 2);
  assert(bond.sites_disjoint());
  assert((bond.type() == "ISING") || (bond.type() == "TJISING"));

  coeff_t J = bond.coupling<coeff_t>();
  int s1 = bond[0];
  int s2 = bond[1];
  bit_t mask = ((bit_t)1 << s1) | ((bit_t)1 << s2);
  bit_t s1_mask = (bit_t)1 << s1;
  bit_t s2_mask = (bit_t)1 << s2;

  // Set values for same/diff (tJ block definition)
  coeff_t val_same, val_diff;
  if (bond.type() == "ISING") {
    val_same = J / 4.;
    val_diff = -J / 4.;
  } else { // (bond.type() == "TJISING")
    val_same = 0.;
    val_diff = -J / 2.;
  }

  auto term_action = [&](bit_t up, bit_t dn) {
    bool b1_up = bitops::gbit(up, s1);
    bool b2_up = bitops::gbit(up, s2);
    bool b1_dn = bitops::gbit(dn, s1);
    bool b2_dn = bitops::gbit(dn, s2);

    if ((b1_up && b2_up) || (b1_dn && b2_dn)) {
      return val_same;
    } else if ((b1_up && b2_dn) || (b1_dn && b2_up)){
      return val_diff;
    } else {
      return 0.;
    }
  };

  tj::generic_term_diag(indexing, term_action, fill);
}

} // namespace hydra::electron
