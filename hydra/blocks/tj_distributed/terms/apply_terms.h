#pragma once

#include <hydra/common.h>

#include <hydra/blocks/tj_distributed/terms/apply_exchange.h>
#include <hydra/blocks/tj_distributed/terms/apply_hopping.h>
#include <hydra/blocks/tj_distributed/terms/apply_ising.h>
#include <hydra/blocks/tj_distributed/terms/apply_number.h>
#include <hydra/blocks/tj_distributed/terms/apply_raise_lower.h>

namespace hydra::tj_distributed {

template <typename bit_t, typename coeff_t, class BasisIn, class BasisOut>
void apply_terms(BondList const &bonds, BasisIn const &basis_in,
                 arma::Col<coeff_t> const &vec_in, BasisOut const &basis_out,
                 arma::Col<coeff_t> &vec_out) try {

  (void)basis_out;
  for (auto bond : bonds) {
    std::string type = bond.type();
    if ((type == "ISING") || (type == "TJISING")) {
      tj_distributed::apply_ising<bit_t, coeff_t>(bond, basis_in, vec_in,
                                                  vec_out);
    } else if (type == "EXCHANGE") {
      tj_distributed::apply_exchange<bit_t, coeff_t>(bond, basis_in, vec_in,
                                                     vec_out);

      // } else if ((type == "NUMBERUP") || (type == "NUMBERDN")) {
      //   tj_distributed::apply_number<bit_t, coeff_t>(bond, basis_in, fill);
      // }
      // else if ((type == "HOPUP") || (type == "HOPDN")) {
      //   tj_distributed::apply_hopping<bit_t, coeff_t>(bond, basis_in, fill);
      // } else if ((type == "CDAGUP") || (type == "CUP") || (type == "CDAGDN")
      // ||
      //            (type == "CDN")) {
      //   tj_distributed::apply_raise_lower<bit_t, coeff_t>(bond, basis_in,
      //                                                     basis_out, fill);
      // }
    } else {
      HydraThrow(std::runtime_error,
                 "Unknown bond type for \"tJDistributed\" block");
    }
  }
} catch (...) {
  HydraRethrow("Unable to apply terms for \"tJDistributed\" block");
}

} // namespace hydra::tj_distributed
