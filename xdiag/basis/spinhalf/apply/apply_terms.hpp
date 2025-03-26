#pragma once

#include <xdiag/basis/apply_identity.hpp>
#include <xdiag/basis/spinhalf/apply/apply_exchange.hpp>
#include <xdiag/basis/spinhalf/apply/apply_matrix.hpp>
#include <xdiag/basis/spinhalf/apply/apply_scalar_chirality.hpp>
#include <xdiag/basis/spinhalf/apply/apply_spsm.hpp>
#include <xdiag/basis/spinhalf/apply/apply_sz.hpp>
#include <xdiag/basis/spinhalf/apply/apply_szsz.hpp>

#include <xdiag/common.hpp>
#include <xdiag/utils/timing.hpp>

namespace xdiag::basis::spinhalf {

template <typename bit_t, typename coeff_t, bool symmetric, class basis_t,
          class fill_f>
void apply_terms(OpSum const &ops, basis_t const &basis_in,
                 basis_t const &basis_out, fill_f fill) try {
  for (auto const &[cpl, op] : ops.plain()) {
    std::string type = op.type();
    if (type == "Id") {
      apply_identity<coeff_t>(cpl, basis_in, fill);
    } else if (type == "Exchange") {
      spinhalf::apply_exchange<bit_t, coeff_t, symmetric>(cpl, op, basis_in,
                                                          basis_out, fill);
    } else if (type == "SzSz") {
      spinhalf::apply_szsz<bit_t, coeff_t, symmetric>(cpl, op, basis_in,
                                                      basis_out, fill);
    } else if (type == "Sz") {
      spinhalf::apply_sz<bit_t, coeff_t, symmetric>(cpl, op, basis_in,
                                                    basis_out, fill);
    } else if (type == "S+") {
      spinhalf::apply_spsm<bit_t, coeff_t, symmetric>(cpl, op, basis_in,
                                                      basis_out, fill);
    } else if (type == "S-") {
      spinhalf::apply_spsm<bit_t, coeff_t, symmetric>(cpl, op, basis_in,
                                                      basis_out, fill);
    } else if (type == "ScalarChirality") {
      spinhalf::apply_scalar_chirality<bit_t, coeff_t, symmetric>(
          cpl, op, basis_in, basis_out, fill);
    } else if (type == "Matrix") {
      spinhalf::apply_matrix<bit_t, coeff_t, symmetric>(cpl, op, basis_in,
                                                        basis_out, fill);
    } else {
      XDIAG_THROW(
          fmt::format("Unknown Op type for Spinhalf block: \"{}\"", type));
    }
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

} // namespace xdiag::basis::spinhalf
