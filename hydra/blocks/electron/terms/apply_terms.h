#pragma once

#include <hydra/bitops/bitops.h>
#include <hydra/common.h>

#include <hydra/blocks/electron/terms/apply_exchange.h>
#include <hydra/blocks/electron/terms/apply_hopping.h>
#include <hydra/blocks/electron/terms/apply_ising.h>
#include <hydra/blocks/electron/terms/apply_u.h>

#include <hydra/blocks/electron/terms/apply_symmetric_exchange.h>
#include <hydra/blocks/electron/terms/apply_symmetric_hopping.h>
#include <hydra/blocks/electron/terms/apply_symmetric_ising.h>
#include <hydra/blocks/electron/terms/apply_symmetric_u.h>

namespace hydra::electron {

template <typename bit_t, typename coeff_t, bool symmetric, class IndexingIn,
          class IndexingOut, class Fill>
void apply_terms(BondList const &bonds, IndexingIn const &indexing_in,
                 IndexingOut const &indexing_out, Fill &fill) {

  if constexpr (symmetric) {

    for (auto bond : bonds) {

      if (bond.type_defined()) {

        if (bond.type() == "EXCHANGE") {
          electron::apply_symmetric_exchange<bit_t, coeff_t>(bond, indexing_in,
                                                             fill);
        } else if (bond.type() == "ISING") {
          electron::apply_symmetric_ising<bit_t, coeff_t>(bond, indexing_in,
                                                          fill);
        } else if (bond.type() == "HOPUP") {
          electron::apply_symmetric_hopping<bit_t, coeff_t>(bond, indexing_in,
                                                            fill);
        } else if (bond.type() == "HOPDN") {
          electron::apply_symmetric_hopping<bit_t, coeff_t>(bond, indexing_in,
                                                            fill);
        } else {
          Log.err("Error in electron::apply_terms: Unknown bond type {}",
                  bond.type());
        }
      }
    }

    if (bonds.coupling_defined("U")) {
      coeff_t U = bonds.coupling<coeff_t>("U");
      electron::apply_symmetric_u<bit_t, coeff_t>(U, indexing_in, fill);
    }

  } else { // not symmetric

    for (auto bond : bonds) {

      if (bond.type_defined()) {

        if (bond.type() == "EXCHANGE") {
          electron::apply_exchange<bit_t, coeff_t>(bond, indexing_in, fill);
        } else if (bond.type() == "ISING") {
          electron::apply_ising<bit_t, coeff_t>(bond, indexing_in, fill);
        } else if (bond.type() == "HOPUP") {
          electron::apply_hopping<bit_t, coeff_t>(bond, indexing_in, fill);
        } else if (bond.type() == "HOPDN") {
          electron::apply_hopping<bit_t, coeff_t>(bond, indexing_in, fill);
        } else {
          Log.err("Error in electron::apply_terms: Unknown bond type {}",
                  bond.type());
        }
      }
    }

    if (bonds.coupling_defined("U")) {
      coeff_t U = bonds.coupling<coeff_t>("U");
      electron::apply_u<bit_t, coeff_t>(U, indexing_in, fill);
    }
  }
}

} // namespace hydra::electron
