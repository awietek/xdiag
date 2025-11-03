// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/basis/tj/apply/apply_exchange.hpp>
#include <xdiag/basis/tj/apply/apply_hopping.hpp>
#include <xdiag/basis/tj/apply/apply_number.hpp>
#include <xdiag/basis/tj/apply/apply_number_number.hpp>
#include <xdiag/basis/tj/apply/apply_raise_lower.hpp>
#include <xdiag/basis/tj/apply/apply_szsz.hpp>
#include <xdiag/common.hpp>

namespace xdiag::basis::tj {

template <typename coeff_t, bool symmetric, class basis_t, class fill_f>
void apply_terms(OpSum const &ops, basis_t const &basis_in,
                 basis_t const &basis_out, fill_f fill) try {

  for (auto const &[cpl, op] : ops) {
    std::string type = op.type();
    if ((type == "SzSz") || (type == "tJSzSz")) {
      tj::apply_szsz<symmetric, coeff_t>(cpl, op, basis_in, basis_out, fill);
    } else if ((type == "Nup") || (type == "Ndn")) {
      tj::apply_number<symmetric, coeff_t>(cpl, op, basis_in, basis_out, fill);
    } else if (type == "NtotNtot") {
      tj::apply_number_number<symmetric, coeff_t>(cpl, op, basis_in, basis_out,
                                                  fill);
    } else if (type == "Exchange") {
      tj::apply_exchange<symmetric, coeff_t>(cpl, op, basis_in, basis_out,
                                             fill);
    } else if ((type == "Hopup") || (type == "Hopdn")) {
      tj::apply_hopping<symmetric, coeff_t>(cpl, op, basis_in, basis_out, fill);
    } else if ((type == "Cdagup") || (type == "Cup") || (type == "Cdagdn") ||
               (type == "Cdn")) {
      tj::apply_raise_lower<symmetric, coeff_t>(cpl, op, basis_in, basis_out,
                                                fill);
    } else {
      XDIAG_THROW(fmt::format("Unknown Op type for tJ block: \"{}\"", type));
    }
  }
}
XDIAG_CATCH

} // namespace xdiag::basis::tj
