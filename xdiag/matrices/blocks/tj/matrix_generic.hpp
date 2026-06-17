// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <string>

#include <xdiag/algebra/algebras/tj_implementation_algebra.hpp>
#include <xdiag/algebra/normal_order.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/operators/valid.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/format.hpp>

#include <xdiag/matrices/blocks/tj/terms/term_cdagc_string.hpp>
#include <xdiag/matrices/blocks/tj/terms/term_exchange.hpp>
#include <xdiag/matrices/blocks/tj/terms/term_hopdn.hpp>
#include <xdiag/matrices/blocks/tj/terms/term_hopup.hpp>
#include <xdiag/matrices/blocks/tj/terms/term_number.hpp>
#include <xdiag/matrices/blocks/tj/terms/term_raise_lower.hpp>
#include <xdiag/matrices/blocks/tj/terms/term_szsz.hpp>
#include <xdiag/matrices/terms/term_identity.hpp>

namespace xdiag::matrices::tj {

template <typename coeff_t, typename basis_t, typename fill_f>
void matrix_generic(OpSum const &ops, basis_t const &basis_in,
                    basis_t const &basis_out, fill_f fill) try {

  operators::check_valid(ops);
  auto algebra = algebra::tj_implementation_algebra(basis_in.nsites());
  auto ops_compiled = normal_order(ops.plain(), algebra);

  for (auto const &[c, monomial] : ops_compiled) {
    // Products of elementary operators are handled by the general Cdag/C string
    // kernel; size-1 named operators dispatch to their dedicated fast kernels.
    if (monomial.size() != 1) {
      term_cdagc_string<coeff_t>(c, monomial, basis_in, basis_out, fill);
      continue;
    }
    Op op = monomial[0];
    std::string type = op.type();

    if (type == "Id") {
      matrices::term_identity<coeff_t>(c, basis_in, fill);
    } else if (type == "Hopup") {
      term_hopup<coeff_t>(c, op, basis_in, basis_out, fill);
    } else if (type == "Hopdn") {
      term_hopdn<coeff_t>(c, op, basis_in, basis_out, fill);
    } else if (type == "Cup") {
      term_cup<coeff_t>(c, op, basis_in, basis_out, fill);
    } else if (type == "Cdagup") {
      term_cdagup<coeff_t>(c, op, basis_in, basis_out, fill);
    } else if (type == "Cdn") {
      term_cdn<coeff_t>(c, op, basis_in, basis_out, fill);
    } else if (type == "Cdagdn") {
      term_cdagdn<coeff_t>(c, op, basis_in, basis_out, fill);
    } else if ((type == "Nup") || (type == "Ndn")) {
      term_number<coeff_t>(c, op, basis_in, basis_out, fill);
    } else if ((type == "NupNup") || (type == "NdnNdn") || (type == "NupNdn") ||
               (type == "NdnNup") || (type == "NtotNtot")) {
      term_number_number<coeff_t>(c, op, basis_in, basis_out, fill);
    } else if ((type == "SzSz") || (type == "tJSzSz")) {
      term_szsz<coeff_t>(c, op, basis_in, basis_out, fill);
    } else if (type == "Exchange") {
      term_exchange<coeff_t>(c, op, basis_in, basis_out, fill);
    } else {
      XDIAG_THROW(fmt::format("Unknown Op type for tJ basis \"{}\"", type));
    }
  }
}
XDIAG_CATCH

} // namespace xdiag::matrices::tj
