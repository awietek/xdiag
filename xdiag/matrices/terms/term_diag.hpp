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

namespace xdiag::matrices {

// Implementation for:
// BasisOnTheFly
// BasisSymmetric

template <typename basis_t, typename term_coeff_f, typename fill_f>
void term_diag(basis_t const &basis, term_coeff_f term_coeff, fill_f fill) {

  // OpenMP parallel implementation
#ifdef _OPENMP
  int64_t size = basis.size();

#pragma omp parallel
  {
    int num_thread = omp_get_thread_num();
    int nthreads = omp_get_num_threads();
    int64_t idx = num_thread * (size / nthreads);
    auto begin = basis.begin() + idx;
    auto end = (num_thread == nthreads - 1)
                   ? basis.end()
                   : basis.begin() + (num_thread + 1) * (size / nthreads);
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

} // namespace xdiag::matrices
