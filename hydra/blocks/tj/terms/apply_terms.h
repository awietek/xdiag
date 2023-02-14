#pragma once

#include <hydra/bitops/bitops.h>
#include <hydra/common.h>

#include <hydra/blocks/tj/terms/apply_exchange.h>
#include <hydra/blocks/tj/terms/apply_hopping.h>
#include <hydra/blocks/tj/terms/apply_ising.h>
#include <hydra/blocks/tj/terms/apply_number.h>

#include <hydra/blocks/tj/terms/apply_symmetric_exchange.h>
#include <hydra/blocks/tj/terms/apply_symmetric_hopping.h>
#include <hydra/blocks/tj/terms/apply_symmetric_ising.h>

namespace hydra::tj {

template <typename bit_t, typename coeff_t, bool symmetric, class IndexingIn,
          class IndexingOut, class Fill>
void apply_terms(BondList const &bonds, IndexingIn const &indexing_in,
                 IndexingOut const &indexing_out, Fill &fill) {
  (void) indexing_out;
  
  if constexpr (symmetric) {

    for (auto bond : bonds) {

      if (bond.type_defined()) {

        if (bond.type() == "EXCHANGE") {
          tj::apply_symmetric_exchange<bit_t, coeff_t>(bond, indexing_in, fill);
        } else if ((bond.type() == "ISING") || (bond.type() == "TJISING")) {
          tj::apply_symmetric_ising<bit_t, coeff_t>(bond, indexing_in, fill);
        } else if (bond.type() == "HOPUP") {
          tj::apply_symmetric_hopping<bit_t, coeff_t>(bond, indexing_in, fill);
        } else if (bond.type() == "HOPDN") {
          tj::apply_symmetric_hopping<bit_t, coeff_t>(bond, indexing_in, fill);
        } else {
          Log.err("Error in tj::apply_terms: Unknown bond type {}",
                  bond.type());
        }
      }
    }

  } else { // not symmetric

    for (auto bond : bonds) {

      if (bond.type_defined()) {

        if (bond.type() == "EXCHANGE") {
          tj::apply_exchange<bit_t, coeff_t>(bond, indexing_in, fill);
        } else if ((bond.type() == "ISING") || (bond.type() == "TJISING")) {
          tj::apply_ising<bit_t, coeff_t>(bond, indexing_in, fill);
        } else if (bond.type() == "HOPUP") {
          tj::apply_hopping<bit_t, coeff_t>(bond, indexing_in, fill);
        } else if (bond.type() == "HOPDN") {
          tj::apply_hopping<bit_t, coeff_t>(bond, indexing_in, fill);
        } else if ((bond.type() == "NUMBERUP") || (bond.type() == "NUMBERDN")) {
          tj::apply_number<bit_t, coeff_t>(bond, indexing_in, fill);
        } else {
          Log.err("Error in tj::apply_terms: Unknown bond type {}",
                  bond.type());
        }
      }
    }
  }
}

} // namespace hydra::tj
