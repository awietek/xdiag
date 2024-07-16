#pragma once

#include <xdiag/blocks/tj/compile.hpp>
#include <xdiag/blocks/tj/terms/apply_exchange.hpp>
#include <xdiag/blocks/tj/terms/apply_hopping.hpp>
#include <xdiag/blocks/tj/terms/apply_ising.hpp>
#include <xdiag/blocks/tj/terms/apply_number.hpp>
#include <xdiag/blocks/tj/terms/apply_raise_lower.hpp>
#include <xdiag/common.hpp>

namespace xdiag::tj {

template <typename bit_t, typename coeff_t, bool symmetric, class BasisIn,
          class BasisOut, class Fill>
void apply_terms(OpSum const &ops, BasisIn const &basis_in,
                 BasisOut const &basis_out, Fill &fill, double zero_precision) {

  OpSum ops_compiled =
      tj::compile(ops, basis_in.n_sites(), zero_precision);
  for (auto op : ops_compiled) {
    std::string type = op.type();
    if ((type == "ISING") || (type == "TJISING")) {
      tj::apply_ising<bit_t, coeff_t, symmetric>(op, basis_in, fill);
    } else if ((type == "NUMBERUP") || (type == "NUMBERDN")) {
      tj::apply_number<bit_t, coeff_t, symmetric>(op, basis_in, fill);
    } else if (type == "EXCHANGE") {
      tj::apply_exchange<bit_t, coeff_t, symmetric>(op, basis_in, fill);
    } else if ((type == "HOPUP") || (type == "HOPDN")) {
      tj::apply_hopping<bit_t, coeff_t, symmetric>(op, basis_in, fill);
    } else if ((type == "CDAGUP") || (type == "CUP") || (type == "CDAGDN") ||
               (type == "CDN")) {
      tj::apply_raise_lower<bit_t, coeff_t, symmetric>(op, basis_in,
                                                       basis_out, fill);
    }
  }
}

} // namespace xdiag::tj
