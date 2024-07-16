#pragma once

#include <xdiag/blocks/spinhalf/terms/apply_exchange.hpp>
#include <xdiag/blocks/spinhalf/terms/apply_ising.hpp>
#include <xdiag/blocks/spinhalf/terms/apply_non_branching.hpp>
#include <xdiag/blocks/spinhalf/terms/apply_scalar_chirality.hpp>
#include <xdiag/blocks/spinhalf/terms/apply_spsm.hpp>
#include <xdiag/blocks/spinhalf/terms/apply_sz.hpp>
#include <xdiag/common.hpp>
#include <xdiag/utils/print_macro.hpp>
#include <xdiag/utils/timing.hpp>

namespace xdiag::spinhalf {

template <typename bit_t, typename coeff_t, bool symmetric, class BasisIn,
          class BasisOut, class Fill>
void apply_terms(OpSum const &ops, BasisIn const &basis_in,
                 BasisOut const &basis_out, Fill &fill,
                 double zero_precision) try {
  OpSum ops_compiled =
      spinhalf::compile(ops, basis_in.n_sites(), zero_precision);

  for (auto op : ops_compiled) {

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

} // namespace xdiag::spinhalf
