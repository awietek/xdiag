#pragma once

#include <hydra/bitops/bitops.h>
#include <hydra/common.h>

#include <hydra/blocks/spinhalf/terms/apply_term_diag.h>

namespace hydra::terms::spinhalf {

// Ising term: H S^z_i

template <typename bit_t, typename coeff_t, class Indexing, class Fill>
void apply_sz(Bond const &bond, Couplings const &couplings, Indexing &&indexing,
              Fill &&fill) {
  if (!(bond.size() == 1)) {
    Log.err("Error in spinhalf::apply_sz: bond has {} sites, expected 1",
            bond.size());
  }

  int s = bond[0];
  bit_t mask = ((bit_t)1 << s);

  std::string cpl = bond.coupling();
  coeff_t H = utils::get_coupling<coeff_t>(couplings, cpl);
  coeff_t val_up = H / 2.;
  coeff_t val_dn = -H / 2.;

  // Function to flip spin and get coefficient
  auto term_coeff = [&mask, &val_up, &val_dn](bit_t spins) -> coeff_t {
    if (spins & mask) {
      return val_up;
    } else {
      return val_dn;
    }
  };

  if (!lila::close(H, 0.)) {
    apply_term_diag<bit_t, coeff_t>(indexing, term_coeff, fill);
  }
}

} // namespace hydra::terms::spinhalf
