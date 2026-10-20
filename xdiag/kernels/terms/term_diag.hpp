// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>

#include <xdiag/basis/basis_onthefly.hpp>
#include <xdiag/kernels/terms/term_offdiag.hpp>
#include <xdiag/utils/thread_range.hpp>

namespace xdiag::kernels {

// Implementation for:
// BasisOnTheFly
// BasisSymmetric

template <typename basis_t, typename term_coeff_f, typename fill_f>
void term_diag(basis_t const &basis_in, basis_t const &basis_out,
               term_coeff_f term_coeff, fill_f fill) {

  // Input and output basis can differ for different irreps,
  // where the action is still diagonal in the basis
  if (basis_in != basis_out) {
    using bit_t = typename basis_t::bit_t;
    term_offdiag(
        basis_in, basis_out, [](bit_t) { return true; },
        [&term_coeff](bit_t s) { return std::pair{s, term_coeff(s)}; }, fill);
  } else {
    basis_t const &basis = basis_in;

    utils::for_each_thread_range(basis, [&](auto begin, auto end, int64_t idx) {
      for (auto it = begin; it != end; ++it, ++idx) {
        auto coeff = term_coeff(*it);
        fill(idx, idx, coeff);
      }
    });
  }
}

} // namespace xdiag::kernels
