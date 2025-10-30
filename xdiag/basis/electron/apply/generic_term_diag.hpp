// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <tuple>

#include <xdiag/basis/electron/apply/generic_term.hpp>
#include <xdiag/parallel/omp/omp_utils.hpp>

namespace xdiag::basis::electron {

template <typename coeff_t, typename basis_t, typename apply_f, typename fill_f>
void generic_term_diag_sym(basis_t const &basis, apply_f apply, fill_f fill) {
  using bit_t = typename basis_t::bit_t;

#ifdef _OPENMP
#pragma omp parallel
  {
    int num_thread = omp_get_thread_num();
#pragma omp for schedule(runtime)
#endif
    for (int64_t idx_ups = 0; idx_ups < basis.n_rep_ups(); ++idx_ups) {
      bit_t ups = basis.rep_ups(idx_ups);
      auto dnss = basis.dns_for_ups_rep(ups);
      int64_t idx = basis.ups_offset(idx_ups);
      for (bit_t dns : dnss) {
        coeff_t c = apply(ups, dns);
        XDIAG_FILL(idx, idx, c);
        ++idx;
      }
    }
#ifdef _OPENMP
  }
#endif
}

template <typename coeff_t, typename basis_t, typename apply_f, typename fill_f>
void generic_term_diag_no_sym(basis_t const &basis, apply_f apply,
                              fill_f fill) {
  using bit_t = typename basis_t::bit_t;

#ifdef _OPENMP
#pragma omp parallel
  {
    int num_thread = omp_get_thread_num();
    auto ups_and_idces = basis.states_indices_ups_thread();
#else
  auto ups_and_idces = basis.states_indices_ups();
#endif

    for (auto [ups, idx_up] : ups_and_idces) {
      int64_t idx = idx_up * basis.size_dns();
      for (bit_t dns : basis.states_dns()) {
        coeff_t c = apply(ups, dns);
        XDIAG_FILL(idx, idx, c);
        ++idx;
      }
    }
#ifdef _OPENMP
  }
#endif
}

template <bool symmetric, typename coeff_t, typename basis_t, typename apply_f,
          typename fill_f>
void generic_term_diag(basis_t const &basis_in, basis_t const &basis_out,
                       apply_f apply, fill_f fill) {
  if (basis_in == basis_out) {
    if constexpr (symmetric) {
      generic_term_diag_sym<coeff_t>(basis_in, apply, fill);
    } else {
      generic_term_diag_no_sym<coeff_t>(basis_in, apply, fill);
    }
  } else { // basis_in != basis_out
    using bit_t = typename basis_t::bit_t;
    auto apply_offdiag = [&](bit_t ups,
                             bit_t dns) -> std::tuple<coeff_t, bit_t, bit_t> {
      return {apply(ups, dns), ups, dns};
    };

    if constexpr (symmetric) {
      generic_term_sym<coeff_t>(basis_in, basis_out, apply_offdiag, fill);
    } else {
      generic_term_no_sym<coeff_t>(basis_in, basis_out, apply_offdiag, fill);
    }
  }
}

} // namespace xdiag::basis::electron
