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
      electron::apply_hopping<coeff_t, symmetric>(cpl, op, basis_in, basis_out,
                                                  fill);
    } else if ((type == "Cdagup") || (type == "Cup") || (type == "Cdagdn") ||
               (type == "Cdn")) {
      electron::apply_raise_lower<coeff_t, symmetric>(cpl, op, basis_in,
                                                      basis_out, fill);
    } else if (type == "SzSz") {
      electron::apply_szsz<coeff_t, symmetric>(cpl, op, basis_in, basis_out,
                                               fill);
    } else if (type == "Exchange") {
      electron::apply_exchange<coeff_t, symmetric>(cpl, op, basis_in, basis_out,
                                                   fill);
    } else if ((type == "Nup") || (type == "Ndn")) {
      electron::apply_number<coeff_t, symmetric>(cpl, op, basis_in, basis_out,
                                                 fill);
    } else if (type == "Nupdn") {
      electron::apply_nupdn<coeff_t, symmetric>(cpl, op, basis_in, basis_out,
                                                fill);
    } else if (type == "NupdnNupdn") {
      electron::apply_nupdn_nupdn<coeff_t, symmetric>(cpl, op, basis_in,
                                                      basis_out, fill);
    } else if (type == "NtotNtot") {
      electron::apply_ntot_ntot<coeff_t, symmetric>(cpl, op, basis_in,
                                                    basis_out, fill);
    } else if (type == "HubbardU") {
      electron::apply_u<coeff_t, symmetric>(cpl, basis_in, fill);
    } else if (type == "NupNdn") {
      electron::apply_nup_ndn<coeff_t, symmetric>(cpl, op, basis_in, basis_out,
                                                  fill);
    } else if (type == "NupNup") {
      electron::apply_nup_nup<coeff_t, symmetric>(cpl, op, basis_in, basis_out,
                                                  fill);
    } else if (type == "NdnNdn") {
      electron::apply_ndn_ndn<coeff_t, symmetric>(cpl, op, basis_in, basis_out,
                                                  fill);
    } else if (type == "Id") {
      apply_identity<coeff_t, basis_t, fill_f>(cpl, basis_in, fill);
    } else {
      XDIAG_THROW(
          fmt::format("Unknown Op type for Electron block: \"{}\"", type));
    }
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

} // namespace xdiag::basis::electron
