// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/matrices/boson/matrix_generic.hpp>
#include <xdiag/operators/opsum.hpp>

namespace xdiag::matrices::boson {

// Matrix policy for the generic kernels in xdiag/matrices/kernels.hpp.
// Wraps boson::matrix_generic as a static callable so it can be passed
// as a template argument: MatrixPolicy::template call<coeff_t>(ops, bi, bo,
// fill).
struct MatrixPolicy {
  template <typename coeff_t, typename basis_t, typename fill_f>
  static void call(OpSum const &ops, basis_t const &basis_in,
                   basis_t const &basis_out, fill_f &&fill) {
    boson::matrix_generic<coeff_t>(ops, basis_in, basis_out, fill);
  }
};

} // namespace xdiag::matrices::boson
