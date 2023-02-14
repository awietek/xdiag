#pragma once

#include <hydra/bitops/bitops.h>
#include <hydra/common.h>

#include <hydra/blocks/spinhalf/terms/apply_term_diag.h>

namespace hydra::spinhalf {

// Ising term: H S^z_i

template <typename bit_t, typename coeff_t, bool symmetric, class IndexingIn,
          class IndexingOut, class Fill>
void apply_sz(Bond const &bond, IndexingIn &&indexing_in,
              IndexingOut &&indexing_out, Fill &&fill) {
  assert(bond.coupling_defined());
  assert(bond.type_defined() && (bond.type() == "SZ"));
  assert(bond.size() == 1);

  int s = bond[0];
  bit_t mask = ((bit_t)1 << s);

  coeff_t H = bond.coupling<coeff_t>();
  coeff_t val_up = H / 2.;
  coeff_t val_dn = -H / 2.;

  if (indexing_in == indexing_out) {

    auto term_coeff = [&mask, &val_up, &val_dn](bit_t spins) -> coeff_t {
      if (spins & mask) {
        return val_up;
      } else {
        return val_dn;
      }
    };

    spinhalf::apply_term_diag<bit_t, coeff_t>(indexing_in, term_coeff, fill);
  } else {
    auto non_zero_term = [](bit_t spins) -> bool {
      return true;
      (void)spins;
    };
    auto term_action = [&mask, &val_up,
                        &val_dn](bit_t spins) -> std::pair<bit_t, coeff_t> {
      if (spins & mask) {
        return {spins, val_up};
      } else {
        return {spins, val_dn};
      }
    };

    if constexpr (symmetric) {
      spinhalf::apply_term_offdiag_sym<bit_t, coeff_t>(
          indexing_in, indexing_out, non_zero_term, term_action, fill);
    } else {
      spinhalf::apply_term_offdiag_no_sym<bit_t, coeff_t>(
          indexing_in, indexing_out, non_zero_term, term_action, fill);
    }
  }
}

} // namespace hydra::spinhalf
