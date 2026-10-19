// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/operators/opsum.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/format.hpp>

#include <xdiag/kernels/terms/term_identity.hpp>
#include <xdiag/kernels/terms/term_matrix.hpp>

namespace xdiag::kernels::boson {

template <typename coeff_t, typename basis_t, typename fill_f>
void matrix_generic(OpSum const &ops, basis_t const &basis_in,
                    basis_t const &basis_out, fill_f fill) try {

  if (basis_in.d() != basis_out.d()) {
    XDIAG_THROW(
        fmt::format("Incompatible local dimensions: input d={}, output d={}",
                    basis_in.d(), basis_out.d()));
  }

  // ops is compiled (see kernels/compile.hpp): normal-ordered with this
  // block's implementation algebra
  for (auto const &[c, monomial] : ops) {
    assert(monomial.size() == 1); // required for properly compiled ops

    Op op = monomial[0];
    std::string type = op.type();

    if (type == "Id") {
      kernels::term_identity<coeff_t>(c, basis_in, fill);
    } else if (type == "Matrix") {
      kernels::term_matrix<coeff_t>(c, op, basis_in, basis_out, fill);
    } else {
      XDIAG_THROW(fmt::format("Unknown Op type for Boson basis \"{}\"", type));
    }
  }
}
XDIAG_CATCH

} // namespace xdiag::kernels::boson
