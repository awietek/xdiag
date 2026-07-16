// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once
#include <cstdint>
#include <xdiag/operators/coeff.hpp>

#include <xdiag/kernels/fill_functions.hpp>

namespace xdiag::kernels {

template <typename coeff_t, class basis_t, class fill_f>
void term_identity(Coeff const &coeff, basis_t const &basis, fill_f fill) try {
  coeff_t s = coeff.scalar().as<coeff_t>();
#ifdef _OPENMP
#pragma omp parallel
  {
    int num_thread = omp_get_thread_num();
#pragma omp for schedule(runtime)
#endif
    for (int64_t i = 0; i < basis.size(); ++i) {
      XDIAG_FILL(i, i, s);
    }
#ifdef _OPENMP
  }
#endif
}
XDIAG_CATCH

} // namespace xdiag::kernels
