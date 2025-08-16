// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/common.hpp>
#include <xdiag/parallel/omp/omp_utils.hpp>

namespace xdiag::basis::spinhalf {

template <typename coeff_t, class basis_t, class term_coeff_f, class fill_f>
void apply_term_diag(basis_t const &basis, term_coeff_f term_coeff,
                     fill_f fill) {
  using bit_t = typename basis_t::bit_t;

#ifdef _OPENMP
  int64_t size = basis.size();

#pragma omp parallel
  {
    int num_thread = omp_get_thread_num();
#pragma omp for schedule(runtime)
    for (int64_t idx = 0; idx < size; ++idx) {
      bit_t spins = basis.state(idx);
      coeff_t coeff = term_coeff(spins);
      XDIAG_FILL(idx, idx, coeff);
    }
  }

#else
  int64_t idx = 0;
  for (auto spins : basis) {
    coeff_t coeff = term_coeff(spins);
    XDIAG_FILL(idx, idx, coeff);
    ++idx;
  }
#endif
}

} // namespace xdiag::basis::spinhalf
