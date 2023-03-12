#pragma once

#include <hydra/bitops/bitops.h>
#include <hydra/common.h>

#include <hydra/blocks/tj/terms/compile.h>

#include <hydra/blocks/tj/terms/apply_exchange.h>
#include <hydra/blocks/tj/terms/apply_hopping.h>
#include <hydra/blocks/tj/terms/apply_ising.h>
#include <hydra/blocks/tj/terms/apply_number.h>
#include <hydra/blocks/tj/terms/apply_raise_lower.h>

namespace hydra::tj {

template <typename bit_t, typename coeff_t, bool symmetric, class IndexingIn,
          class IndexingOut, class Fill>
void apply_terms(BondList const &bonds, IndexingIn const &indexing_in,
                 IndexingOut const &indexing_out, Fill &fill) {

  BondList bonds_compiled = tj::compile(bonds);

  for (auto bond : bonds_compiled) {
    std::string type = bond.type();
    if ((type == "ISING") || (type == "TJISING")) {
      tj::apply_ising<bit_t, coeff_t, symmetric>(bond, indexing_in, fill);
    } else if ((type == "NUMBERUP") || (type == "NUMBERDN")) {
      tj::apply_number<bit_t, coeff_t, symmetric>(bond, indexing_in, fill);
    } else if (type == "EXCHANGE") {
      tj::apply_exchange<bit_t, coeff_t, symmetric>(bond, indexing_in, fill);
    } else if ((type == "HOPUP") || (type == "HOPDN")) {
      tj::apply_hopping<bit_t, coeff_t, symmetric>(bond, indexing_in, fill);
    } else if ((type == "CDAGUP") || (type == "CUP") || (type == "CDAGDN") ||
               (type == "CDN")) {
      tj::apply_raise_lower<bit_t, coeff_t, symmetric>(bond, indexing_in,
                                                       indexing_out, fill);
    }
  }
}

} // namespace hydra::tj
