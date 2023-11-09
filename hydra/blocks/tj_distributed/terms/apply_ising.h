#pragma once

#include <hydra/common.h>
#include <hydra/operators/bond.h>

#include <hydra/blocks/tj_distributed/terms/generic_term_diag.h>

namespace hydra::tj_distributed {

template <typename bit_t, typename coeff_t, class Basis>
void apply_ising(Bond const &bond, Basis &&basis,
                 const coeff_t* vec_in,
                 coeff_t * vec_out) try {
  assert(bond.coupling_defined());
  assert(bond.type_defined());
  assert(bond.size() == 2);
  assert(bond.sites_disjoint());

  std::string type = bond.type();
  assert((type == "ISING") || (type == "TJISING"));

  coeff_t J = bond.coupling<coeff_t>();
  int64_t s1 = bond[0];
  int64_t s2 = bond[1];
  bit_t s1_mask = (bit_t)1 << s1;
  bit_t s2_mask = (bit_t)1 << s2;

  // Set values for same/diff (tJ block definition)
  coeff_t val_same, val_diff;
  if (type == "ISING") {
    val_same = J / 4.;
    val_diff = -J / 4.;
  } else { // (type == "TJISING")
    val_same = 0.;
    val_diff = -J / 2.;
  }

  auto term_action = [&](bit_t up, bit_t dn) -> coeff_t {
    bool b1_up = up & s1_mask;
    bool b2_up = up & s2_mask;
    bool b1_dn = dn & s1_mask;
    bool b2_dn = dn & s2_mask;

    if ((b1_up && b2_up) || (b1_dn && b2_dn)) {
      return val_same;
    } else if ((b1_up && b2_dn) || (b1_dn && b2_up)) {
      return val_diff;
    } else {
      return 0.;
    }
  };

  tj_distributed::generic_term_diag<bit_t, coeff_t>(basis, term_action, vec_in,
                                                    vec_out);
} catch (...) {
  HydraRethrow("Unable to apply Ising term form \"tJDistributed\" block");
}

} // namespace hydra::tj_distributed
