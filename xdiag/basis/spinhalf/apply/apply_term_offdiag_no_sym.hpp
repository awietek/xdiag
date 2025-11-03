// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/common.hpp>
#include <xdiag/parallel/omp/omp_utils.hpp>

namespace xdiag::basis::spinhalf {

template <typename coeff_t, class basis_t, class non_zero_term_f,
          class term_action_f, class fill_f>
void apply_term_offdiag_no_sym(basis_t const &basis_in,
                               basis_t const &basis_out,
                               non_zero_term_f non_zero_term,
                               term_action_f term_action, fill_f fill) {
  using bit_t = typename basis_t::bit_t;

#ifdef _OPENMP
  int64_t size = basis_in.size();

#pragma omp parallel
  {
    int num_thread = omp_get_thread_num();
#pragma omp for schedule(runtime)
    for (int64_t idx_in = 0; idx_in < size; ++idx_in) {
      bit_t spins_in = basis_in.state(idx_in);
      if (non_zero_term(spins_in)) {
        auto [spins_out, coeff] = term_action(spins_in);
        auto idx_out = basis_out.index(spins_out);
        XDIAG_FILL(idx_in, idx_out, coeff);
      }
    }
  }

#else
  int64_t idx_in = 0;
  for (auto spins_in : basis_in) {
    if (non_zero_term(spins_in)) {
      auto [spins_out, coeff] = term_action(spins_in);
      auto idx_out = basis_out.index(spins_out);
      XDIAG_FILL(idx_in, idx_out, coeff);
    }
    ++idx_in;
  }
#endif
}

} // namespace xdiag::basis::spinhalf
