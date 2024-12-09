#pragma once

#include <xdiag/basis/electron/apply/apply_exchange.hpp>
#include <xdiag/basis/electron/apply/apply_hopping.hpp>
#include <xdiag/basis/electron/apply/apply_ising.hpp>
#include <xdiag/basis/electron/apply/apply_number.hpp>
#include <xdiag/basis/electron/apply/apply_raise_lower.hpp>
#include <xdiag/basis/electron/apply/apply_u.hpp>

#include <xdiag/common.hpp>
#include <xdiag/operators/opsum.hpp>

namespace xdiag::basis::electron {

template <typename bit_t, typename coeff_t, bool symmetric, class BasisIn,
          class BasisOut, class Fill>
void apply_terms(OpSum const &ops, BasisIn const &basis_in,
                 BasisOut const &basis_out, Fill &fill) try {

  for (auto const &[cpl, op] : ops.plain()) {

    std::string type = op.type();
    if ((type == "HOPUP") || (type == "HOPDN")) {
      electron::apply_hopping<bit_t, coeff_t, symmetric>(cpl, op, basis_in,
                                                         fill);
    } else if ((type == "CDAGUP") || (type == "CUP") || (type == "CDAGDN") ||
               (type == "CDN")) {
      electron::apply_raise_lower<bit_t, coeff_t, symmetric>(cpl, op, basis_in,
                                                             basis_out, fill);
    } else if (type == "ISING") {
      electron::apply_ising<bit_t, coeff_t, symmetric>(cpl, op, basis_in, fill);
    } else if ((type == "NUMBERUP") || (type == "NUMBERDN")) {
      electron::apply_number<bit_t, coeff_t, symmetric>(cpl, op, basis_in,
                                                        fill);
    } else if (type == "HUBBARDU") {
      electron::apply_u<bit_t, coeff_t, symmetric>(cpl, basis_in, fill);
    } else {
      XDIAG_THROW(fmt::format("Unknown Op type \"{}\"", type));
    }
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

} // namespace xdiag::basis::electron
