#pragma once

#include <xdiag/basis/spinhalf/apply/apply_exchange.hpp>
#include <xdiag/basis/spinhalf/apply/apply_ising.hpp>
#include <xdiag/basis/spinhalf/apply/apply_non_branching.hpp>
#include <xdiag/basis/spinhalf/apply/apply_scalar_chirality.hpp>
#include <xdiag/basis/spinhalf/apply/apply_spsm.hpp>
#include <xdiag/basis/spinhalf/apply/apply_sz.hpp>
#include <xdiag/common.hpp>
#include <xdiag/utils/timing.hpp>

namespace xdiag::basis::spinhalf {

template <typename bit_t, typename coeff_t, bool symmetric, class BasisIn,
          class BasisOut, class Fill>
void apply_terms(OpSum const &ops, BasisIn const &basis_in,
                 BasisOut const &basis_out, Fill &fill) try {
  for (auto op : ops) {

    if (op.type() == "EXCHANGE") {
      spinhalf::apply_exchange<bit_t, coeff_t, symmetric>(op, basis_in,
                                                          basis_out, fill);
    } else if (op.type() == "ISING") {
      spinhalf::apply_ising<bit_t, coeff_t, symmetric>(op, basis_in, basis_out,
                                                       fill);
    } else if (op.type() == "SZ") {
      spinhalf::apply_sz<bit_t, coeff_t, symmetric>(op, basis_in, basis_out,
                                                    fill);
    } else if (op.type() == "S+") {
      spinhalf::apply_spsm<bit_t, coeff_t, symmetric>(op, basis_in, basis_out,
                                                      fill);
    } else if (op.type() == "S-") {
      spinhalf::apply_spsm<bit_t, coeff_t, symmetric>(op, basis_in, basis_out,
                                                      fill);
    } else if (op.type() == "SCALARCHIRALITY") {
      spinhalf::apply_scalar_chirality<bit_t, coeff_t, symmetric>(
          op, basis_in, basis_out, fill);
    } else if (op.type() == "NONBRANCHINGOP") {
      spinhalf::apply_non_branching<bit_t, coeff_t, symmetric>(op, basis_in,
                                                               basis_out, fill);
    } else {
      XDIAG_THROW(fmt::format(
          "Error in spinhalf::apply_terms: Unknown Op type \"{}\"", op.type()));
    }
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

} // namespace xdiag::basis::spinhalf
