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
void apply_terms(BondList const &bonds, BasisIn const &basis_in,
                 BasisOut const &basis_out, Fill &fill) try {
  for (auto bond : bonds) {

    if (bond.type_defined()) {
      if (bond.type() == "EXCHANGE") {
        spinhalf::apply_exchange<bit_t, coeff_t, symmetric>(bond, basis_in,
                                                            basis_out, fill);
      } else if (bond.type() == "ISING") {
        spinhalf::apply_ising<bit_t, coeff_t, symmetric>(bond, basis_in,
                                                         basis_out, fill);
      } else if (bond.type() == "SZ") {
        spinhalf::apply_sz<bit_t, coeff_t, symmetric>(bond, basis_in, basis_out,
                                                      fill);
      } else if (bond.type() == "S+") {
        spinhalf::apply_spsm<bit_t, coeff_t, symmetric>(bond, basis_in,
                                                        basis_out, fill);
      } else if (bond.type() == "S-") {
        spinhalf::apply_spsm<bit_t, coeff_t, symmetric>(bond, basis_in,
                                                        basis_out, fill);
      } else if (bond.type() == "SCALARCHIRALITY") {
        spinhalf::apply_scalar_chirality<bit_t, coeff_t, symmetric>(
            bond, basis_in, basis_out, fill);
      } else {
        XDIAG_THROW(fmt::format("Error in spinhalf::apply_terms: Unknown bond type \"{}\"", bond.type()));
      }
    } else {
      spinhalf::apply_non_branching<bit_t, coeff_t, symmetric>(bond, basis_in,
                                                               basis_out, fill);
    }
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

} // namespace xdiag::spinhalf
