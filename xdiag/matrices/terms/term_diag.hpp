// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>
#ifdef _OPENMP
#include <omp.h>
#endif

#include <xdiag/basis/basis_onthefly.hpp>
#include <xdiag/matrices/fill_functions.hpp>
#include <xdiag/matrices/terms/term_offdiag.hpp>
#include <xdiag/utils/thread_range.hpp>

namespace xdiag::matrices {

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

    // OpenMP parallel implementation
#ifdef _OPENMP
#pragma omp parallel
    {
      int num_thread = omp_get_thread_num();
      auto [begin, end, idx] =
          utils::thread_range(basis, num_thread, omp_get_num_threads());
      for (auto it = begin; it != end; ++it, ++idx) {
        auto coeff = term_coeff(*it);
        XDIAG_FILL(idx, idx, coeff);
      }
    }

    // Serial implementation
#else
    int64_t idx = 0;
    for (auto spins : basis) {
      auto coeff = term_coeff(spins);
      fill(idx, idx, coeff);
      ++idx;
    }
#endif
  }
}

} // namespace xdiag::matrices
