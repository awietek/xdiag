// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/algebra/algebra.hpp>
#include <xdiag/algebra/normal_order.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/operators/valid.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/format.hpp>
#include <xdiag/utils/xdiag_show.hpp>

#include <xdiag/matrices/blocks/fermion/terms/term_c.hpp>
#include <xdiag/matrices/blocks/fermion/terms/term_cdag.hpp>
#include <xdiag/matrices/blocks/fermion/terms/term_cdagc_string.hpp>
#include <xdiag/matrices/blocks/fermion/terms/term_hop.hpp>
#include <xdiag/matrices/blocks/fermion/terms/term_hop_asym.hpp>
#include <xdiag/matrices/blocks/fermion/terms/term_n.hpp>
#include <xdiag/matrices/blocks/fermion/terms/term_nn.hpp>
#include <xdiag/matrices/blocks/fermion/terms/term_totaln.hpp>
#include <xdiag/matrices/terms/term_identity.hpp>

namespace xdiag::matrices::fermion {

template <typename coeff_t, typename basis_t, typename fill_f>
void matrix_generic(OpSum const &ops, basis_t const &basis_in,
                    basis_t const &basis_out, fill_f fill) try {

  // Get OpSum into format that can be processed
  operators::check_valid(ops);
  auto algebra = algebra::fermion_implementation_algebra(basis_in.nsites());
  auto ops_compiled = normal_order(ops.plain(), algebra);

  for (auto const &[c, monomial] : ops_compiled) {
    if (monomial.size() == 1) {
      Op op = monomial[0];
      std::string type = op.type();

      if (type == "Id") {
        matrices::term_identity<coeff_t>(c, basis_in, fill);
      } else if (type == "Hop") {
        fermion::term_hop<coeff_t>(c, op, basis_in, basis_out, fill);
      } else if (type == "HopAsym") {
        fermion::term_hop_asym<coeff_t>(c, op, basis_in, basis_out, fill);
      } else if (type == "NN") {
        fermion::term_nn<coeff_t>(c, op, basis_in, basis_out, fill);
      } else if (type == "N") {
        fermion::term_n<coeff_t>(c, op, basis_in, basis_out, fill);
      } else if (type == "TotalN") {
        fermion::term_totaln<coeff_t>(c, basis_in, basis_out, fill);
      } else if (type == "Cdag") {
        fermion::term_cdag<coeff_t>(c, op, basis_in, basis_out, fill);
      } else if (type == "C") {
        fermion::term_c<coeff_t>(c, op, basis_in, basis_out, fill);
      } else {
        XDIAG_THROW(
            fmt::format("Unknown Op type for Fermion basis \"{}\"", type));
      }
    } else {
      fermion::term_cdagc_string<coeff_t>(c, monomial, basis_in, basis_out,
                                          fill);
    }
  }
}
XDIAG_CATCH

} // namespace xdiag::matrices::fermion
