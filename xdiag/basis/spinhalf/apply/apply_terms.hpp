#pragma once

#include <xdiag/basis/spinhalf/apply/apply_exchange.hpp>
// #include <xdiag/basis/spinhalf/apply/apply_non_branching.hpp>
#include <xdiag/basis/spinhalf/apply/apply_scalar_chirality.hpp>
#include <xdiag/basis/spinhalf/apply/apply_spsm.hpp>
#include <xdiag/basis/spinhalf/apply/apply_sz.hpp>
#include <xdiag/basis/spinhalf/apply/apply_szsz.hpp>
#include <xdiag/common.hpp>
#include <xdiag/utils/timing.hpp>

namespace xdiag::basis::spinhalf {

template <typename bit_t, typename coeff_t, bool symmetric, class BasisIn,
          class BasisOut, class Fill>
void apply_terms(OpSum const &ops, BasisIn const &basis_in,
                 BasisOut const &basis_out, Fill &fill) try {
  for (auto const &[cpl, op] : ops.plain()) {

    std::string type = op.type();
    if (type == "EXCHANGE") {
      spinhalf::apply_exchange<bit_t, coeff_t, symmetric>(cpl, op, basis_in,
                                                          basis_out, fill);
    } else if (type == "SZSZ") {
      spinhalf::apply_szsz<bit_t, coeff_t, symmetric>(cpl, op, basis_in,
                                                      basis_out, fill);
    } else if (type == "SZ") {
      spinhalf::apply_sz<bit_t, coeff_t, symmetric>(cpl, op, basis_in,
                                                    basis_out, fill);
    } else if (type == "S+") {
      spinhalf::apply_spsm<bit_t, coeff_t, symmetric>(cpl, op, basis_in,
                                                      basis_out, fill);
    } else if (type == "S-") {
      spinhalf::apply_spsm<bit_t, coeff_t, symmetric>(cpl, op, basis_in,
                                                      basis_out, fill);
    } else if (type == "SCALARCHIRALITY") {
      spinhalf::apply_scalar_chirality<bit_t, coeff_t, symmetric>(
          cpl, op, basis_in, basis_out, fill);
    } //  else if (type == "NONBRANCHINGOP") {
    //   spinhalf::apply_non_branching<bit_t, coeff_t, symmetric>(
    //       cpl, op, basis_in, basis_out, fill);
    // }
    else {
      XDIAG_THROW(fmt::format("Unknown Op type \"{}\"", type));
    }
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

} // namespace xdiag::basis::spinhalf
