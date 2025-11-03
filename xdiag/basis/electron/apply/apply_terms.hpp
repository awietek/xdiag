// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/common.hpp>
#include <xdiag/operators/opsum.hpp>

#include <xdiag/basis/apply_identity.hpp>
#include <xdiag/basis/electron/apply/apply_exchange.hpp>
#include <xdiag/basis/electron/apply/apply_hopping.hpp>
#include <xdiag/basis/electron/apply/apply_number.hpp>
#include <xdiag/basis/electron/apply/apply_number_number.hpp>
#include <xdiag/basis/electron/apply/apply_raise_lower.hpp>
#include <xdiag/basis/electron/apply/apply_szsz.hpp>
#include <xdiag/basis/electron/apply/apply_u.hpp>

namespace xdiag::basis::electron {

template <typename coeff_t, bool symmetric, class basis_t, class fill_f>
void apply_terms(OpSum const &ops, basis_t const &basis_in,
                 basis_t const &basis_out, fill_f fill) try {

  for (auto const &[cpl, op] : ops.plain()) {

    std::string type = op.type();
    if ((type == "Hopup") || (type == "Hopdn")) {
      apply_hopping<symmetric, coeff_t>(cpl, op, basis_in, basis_out, fill);
    } else if ((type == "Cdagup") || (type == "Cup") || (type == "Cdagdn") ||
               (type == "Cdn")) {
      apply_raise_lower<symmetric, coeff_t>(cpl, op, basis_in, basis_out, fill);
    } else if (type == "SzSz") {
      apply_szsz<symmetric, coeff_t>(cpl, op, basis_in, basis_out, fill);
    } else if (type == "Exchange") {
      apply_exchange<symmetric, coeff_t>(cpl, op, basis_in, basis_out, fill);
    } else if (type == "Nup") {
      apply_nup<symmetric, coeff_t>(cpl, op, basis_in, basis_out, fill);
    } else if (type == "Ndn") {
      apply_ndn<symmetric, coeff_t>(cpl, op, basis_in, basis_out, fill);
    } else if (type == "Nupdn") {
      apply_nupdn<symmetric, coeff_t>(cpl, op, basis_in, basis_out, fill);
    } else if (type == "NupdnNupdn") {
      apply_nupdn_nupdn<symmetric, coeff_t>(cpl, op, basis_in, basis_out, fill);
    } else if (type == "NtotNtot") {
      apply_ntot_ntot<symmetric, coeff_t>(cpl, op, basis_in, basis_out, fill);
    } else if (type == "HubbardU") {
      apply_u<symmetric, coeff_t>(cpl, basis_in, fill);
    } else if (type == "NupNdn") {
      apply_nup_ndn<symmetric, coeff_t>(cpl, op, basis_in, basis_out, fill);
    } else if (type == "NupNup") {
      apply_nup_nup<symmetric, coeff_t>(cpl, op, basis_in, basis_out, fill);
    } else if (type == "NdnNdn") {
      apply_ndn_ndn<symmetric, coeff_t>(cpl, op, basis_in, basis_out, fill);
    } else if (type == "Id") {
      apply_identity<coeff_t, basis_t, fill_f>(cpl, basis_in, fill);
    } else {
      XDIAG_THROW(
          fmt::format("Unknown Op type for Electron block: \"{}\"", type));
    }
  }
}
XDIAG_CATCH

} // namespace xdiag::basis::electron
