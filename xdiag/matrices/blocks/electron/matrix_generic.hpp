// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <string>

#include <xdiag/algebra/algebra.hpp>
#include <xdiag/algebra/normal_order.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/operators/valid.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/format.hpp>

#include <xdiag/matrices/blocks/electron/terms/term_cdagc_string.hpp>
#include <xdiag/matrices/blocks/electron/terms/term_hop.hpp>
#include <xdiag/matrices/blocks/electron/terms/term_hop_asym.hpp>
#include <xdiag/matrices/blocks/electron/terms/term_number.hpp>
#include <xdiag/matrices/blocks/electron/terms/term_number_number.hpp>
#include <xdiag/matrices/blocks/electron/terms/term_raise_lower.hpp>
#include <xdiag/matrices/blocks/electron/terms/term_szsz.hpp>
#include <xdiag/matrices/blocks/electron/terms/term_u.hpp>
#include <xdiag/matrices/terms/term_identity.hpp>

namespace xdiag::matrices::electron {

template <typename coeff_t, typename basis_t, typename fill_f>
void matrix_generic(OpSum const &ops, basis_t const &basis_in,
                    basis_t const &basis_out, fill_f fill) try {

  operators::check_valid(ops);
  auto algebra = algebra::electron_implementation_algebra(basis_in.nsites());
  auto ops_compiled = normal_order(ops.plain(), algebra);

  for (auto const &[c, monomial] : ops_compiled) {
    // Products of elementary operators (e.g. S+ = Cdagup Cdn) are handled by the
    // Cdag/C string kernel; size-1 named operators dispatch to their kernels.
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
    } else if (type == "HopupAsym") {
      term_hopup_asym<coeff_t>(c, op, basis_in, basis_out, fill);
    } else if (type == "HopdnAsym") {
      term_hopdn_asym<coeff_t>(c, op, basis_in, basis_out, fill);
    } else if (type == "Cup") {
      term_cup<coeff_t>(c, op, basis_in, basis_out, fill);
    } else if (type == "Cdagup") {
      term_cdagup<coeff_t>(c, op, basis_in, basis_out, fill);
    } else if (type == "Cdn") {
      term_cdn<coeff_t>(c, op, basis_in, basis_out, fill);
    } else if (type == "Cdagdn") {
      term_cdagdn<coeff_t>(c, op, basis_in, basis_out, fill);
    } else if (type == "Nup") {
      term_nup<coeff_t>(c, op, basis_in, basis_out, fill);
    } else if (type == "Ndn") {
      term_ndn<coeff_t>(c, op, basis_in, basis_out, fill);
    } else if (type == "Nupdn") {
      term_nupdn<coeff_t>(c, op, basis_in, basis_out, fill);
    } else if (type == "NupNup") {
      term_nupnup<coeff_t>(c, op, basis_in, basis_out, fill);
    } else if (type == "NdnNdn") {
      term_ndnndn<coeff_t>(c, op, basis_in, basis_out, fill);
    } else if (type == "NupNdn") {
      term_nupndn<coeff_t>(c, op, basis_in, basis_out, fill);
    } else if (type == "NdnNup") {
      term_ndnnup<coeff_t>(c, op, basis_in, basis_out, fill);
    } else if (type == "NtotNtot") {
      term_ntotntot<coeff_t>(c, op, basis_in, basis_out, fill);
    } else if (type == "NupdnNupdn") {
      term_nupdnnupdn<coeff_t>(c, op, basis_in, basis_out, fill);
    } else if (type == "SzSz") {
      term_szsz<coeff_t>(c, op, basis_in, basis_out, fill);
    } else if (type == "HubbardU") {
      term_hubbardu<coeff_t>(c, basis_in, basis_out, fill);
    } else {
      XDIAG_THROW(
          fmt::format("Unknown Op type for Electron basis \"{}\"", type));
    }
  }
}
XDIAG_CATCH

} // namespace xdiag::matrices::electron
