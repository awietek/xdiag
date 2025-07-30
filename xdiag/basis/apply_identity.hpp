// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once
#include <xdiag/common.hpp>
#include <xdiag/operators/coupling.hpp>
#include <xdiag/parallel/omp/omp_utils.hpp>

namespace xdiag::basis {

template <typename coeff_t, class basis_t, class fill_f>
void apply_identity(Coupling const &cpl, basis_t const &basis,
                    fill_f fill) try {
  coeff_t s = cpl.scalar().as<coeff_t>();
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
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <typename coeff_t, class basis_t>
void apply_identity_distributed(Coupling const &cpl, const coeff_t *vec_in,
                                coeff_t *vec_out, int64_t size) {
  using bit_t = typename basis_t::bit_t;
  coeff_t s = cpl.scalar().as<coeff_t>();
  for (int64_t idx = 0; idx < size; ++idx) {
    vec_out[idx] += s * vec_in[idx];
  }
}

} // namespace xdiag::basis
