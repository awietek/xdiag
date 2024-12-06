#pragma once

#include <xdiag/basis/tj/apply/apply_exchange.hpp>
#include <xdiag/basis/tj/apply/apply_hopping.hpp>
#include <xdiag/basis/tj/apply/apply_ising.hpp>
#include <xdiag/basis/tj/apply/apply_number.hpp>
#include <xdiag/basis/tj/apply/apply_raise_lower.hpp>
#include <xdiag/common.hpp>

namespace xdiag::basis::tj {

template <typename bit_t, typename coeff_t, bool symmetric, class BasisIn,
          class BasisOut, class Fill>
void apply_terms(OpSum const &ops, BasisIn const &basis_in,
                 BasisOut const &basis_out, Fill &fill) try {

  for (auto const &[cpl, op] : ops) {
    std::string type = op.type();
    if ((type == "ISING") || (type == "TJISING")) {
      apply_ising<bit_t, coeff_t, symmetric>(cpl, op, basis_in, fill);
    } else if ((type == "NUMBERUP") || (type == "NUMBERDN")) {
      apply_number<bit_t, coeff_t, symmetric>(cpl, op, basis_in, fill);
    } else if (type == "EXCHANGE") {
      apply_exchange<bit_t, coeff_t, symmetric>(cpl, op, basis_in, fill);
    } else if ((type == "HOPUP") || (type == "HOPDN")) {
      apply_hopping<bit_t, coeff_t, symmetric>(cpl, op, basis_in, fill);
    } else if ((type == "CDAGUP") || (type == "CUP") || (type == "CDAGDN") ||
               (type == "CDN")) {
      apply_raise_lower<bit_t, coeff_t, symmetric>(cpl, op, basis_in, basis_out,
                                                   fill);
    }
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

} // namespace xdiag::basis::tj
