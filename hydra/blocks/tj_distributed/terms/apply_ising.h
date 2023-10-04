#pragma once
#ifdef HYDRA_ENABLE_MPI

#include <hydra/common.h>
#include <hydra/operators/bond.h>

#include <hydra/blocks/tj/terms/generic_term_diag.h>

namespace hydra::tj_distributed {

template <typename bit_t, typename coeff_t, bool symmetric, class Basis,
          class Fill>
void apply_ising(Bond const &bond, Basis &&basis, Fill &&fill) {
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

  
}

} // namespace hydra::tj_distributed
#endif