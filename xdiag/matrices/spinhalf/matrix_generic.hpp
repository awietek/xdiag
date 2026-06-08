// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/algebra/algebras/spinhalf_implementation_algebra.hpp>
#include <xdiag/algebra/normal_order.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/operators/valid.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/format.hpp>

#include <xdiag/matrices/spinhalf/terms/term_exchange.hpp>
#include <xdiag/matrices/spinhalf/terms/term_exchange_asym.hpp>
#include <xdiag/matrices/spinhalf/terms/term_scalar_chirality.hpp>
#include <xdiag/matrices/spinhalf/terms/term_spsm.hpp>
#include <xdiag/matrices/spinhalf/terms/term_sz.hpp>
#include <xdiag/matrices/spinhalf/terms/term_szsz.hpp>
#include <xdiag/matrices/terms/term_identity.hpp>
#include <xdiag/matrices/terms/term_matrix.hpp>

namespace xdiag::matrices::spinhalf {

template <typename coeff_t, typename basis_t, typename fill_f>
void matrix_generic(OpSum const &ops, basis_t const &basis_in,
                    basis_t const &basis_out, fill_f fill) try {

  // Get OpSum into format that can be processed
  operators::check_valid(ops);
  auto algebra = algebra::spinhalf_implementation_algebra();
  auto ops_compiled = normal_order(ops.plain(), algebra);

  for (auto const &[c, monomial] : ops_compiled) {
    assert(monomial.size() == 1); // required for properly compiled ops

    Op op = monomial[0];
    std::string type = op.type();

    if (type == "Id") {
      matrices::term_identity<coeff_t>(c, basis_in, fill);
    } else if (type == "Exchange") {
      spinhalf::term_exchange<coeff_t>(c, op, basis_in, basis_out, fill);
    } else if (type == "ExchangeAsym") {
      spinhalf::term_exchange_asym<coeff_t>(c, op, basis_in, basis_out, fill);
    } else if (type == "SzSz") {
      spinhalf::term_szsz<coeff_t>(c, op, basis_in, basis_out, fill);
    } else if (type == "Sz") {
      spinhalf::term_sz<coeff_t>(c, op, basis_in, basis_out, fill);
    } else if (type == "S+") {
      spinhalf::term_spsm<coeff_t>(c, op, basis_in, basis_out, fill);
    } else if (type == "S-") {
      spinhalf::term_spsm<coeff_t>(c, op, basis_in, basis_out, fill);
    } else if (type == "ScalarChirality") {
      spinhalf::term_scalar_chirality<coeff_t>(c, op, basis_in, basis_out,
                                               fill);
    } else if (type == "Matrix") {
      matrices::term_matrix<coeff_t>(c, op, basis_in, basis_out, fill);
    } else {
      XDIAG_THROW(
          fmt::format("Unknown Op type for Spinhalf basis \"{}\"", type));
    }
  }
}
XDIAG_CATCH

} // namespace xdiag::matrices::spinhalf
