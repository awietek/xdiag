// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/matrices/blocks/spinhalf/matrix_generic.hpp>
#include <xdiag/operators/opsum.hpp>

namespace xdiag::matrices::spinhalf {

// Matrix policy for the generic kernels in xdiag/matrices/kernels.hpp.
// Wraps spinhalf::matrix_generic as a static callable so it can be passed
// as a template argument: MatrixPolicy::template call<coeff_t>(ops, bi, bo,
// fill).
struct MatrixPolicy {
  template <typename coeff_t, typename basis_t, typename fill_f>
  static void call(OpSum const &ops, basis_t const &basis_in,
                   basis_t const &basis_out, fill_f &&fill) {
    spinhalf::matrix_generic<coeff_t>(ops, basis_in, basis_out, fill);
  }
};

} // namespace xdiag::matrices::spinhalf
