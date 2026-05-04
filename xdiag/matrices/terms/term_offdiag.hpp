// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>
#include <type_traits>
#ifdef _OPENMP
#include <omp.h>
#endif

#include <xdiag/basis/basis_onthefly.hpp>
#include <xdiag/matrices/fill_functions.hpp>
#include <xdiag/utils/likely.hpp>

namespace xdiag::matrices {

template <typename enumeration_t, typename non_zero_term_f,
          typename term_action_f, typename fill_f>
void term_offdiag(basis::BasisOnTheFly<enumeration_t> const &basis_in,
                  basis::BasisOnTheFly<enumeration_t> const &basis_out,
                  non_zero_term_f non_zero_term, term_action_f term_action,
                  fill_f fill) {
  using bit_t = typename enumeration_t::bit_t;

  // OpenMP parallel implementation
#ifdef _OPENMP
#pragma omp parallel
  {
    int num_thread = omp_get_thread_num();
    int nthreads = omp_get_num_threads();
    int64_t size = basis_in.size();
    int64_t idx_in = num_thread * (size / nthreads);
    auto begin = basis_in.begin() + idx_in;
    auto end = (num_thread == nthreads - 1)
                   ? basis_in.end()
                   : basis_in.begin() + (num_thread + 1) * (size / nthreads);
    for (auto it = begin; it != end; ++it, ++idx_in) {
      bit_t spins_in = *it;
      if (non_zero_term(spins_in)) {
        auto [spins_out, coeff] = term_action(spins_in);
        int64_t idx_out = basis_out.index(spins_out);
        XDIAG_FILL(idx_in, idx_out, coeff);
      }
    }
  }

  // Serial implementation
#else
  int64_t idx_in = 0;
  for (auto const &spins_in : basis_in) {
    if (non_zero_term(spins_in)) {
      auto [spins_out, coeff] = term_action(spins_in);
      int64_t idx_out = basis_out.index(spins_out);
      fill(idx_in, idx_out, coeff);
    }
    ++idx_in;
  }
#endif
}

template <typename enumeration_t, typename non_zero_term_f,
          typename term_action_f, typename fill_f>
void term_offdiag(basis::BasisSymmetric<enumeration_t> const &basis_in,
                  basis::BasisSymmetric<enumeration_t> const &basis_out,
                  non_zero_term_f non_zero_term, term_action_f term_action,
                  fill_f fill) {
  using bit_t = typename enumeration_t::bit_t;
  using coeff_t =
      typename std::invoke_result_t<decltype(term_action), bit_t>::second_type;

  Representation irrep = basis_out.irrep();
  arma::Col<coeff_t> characters = irrep.characters().as<arma::Col<coeff_t>>();

#ifdef _OPENMP
#pragma omp parallel
  {
    int num_thread = omp_get_thread_num(); // needed for XDIAG_FILL in case of
                                           // sparse matrices
#pragma omp for schedule(runtime)
    for (int64_t idx_in = 0; idx_in < basis_in.size(); ++idx_in) {
      bit_t spins_in = basis_in[idx_in];
      if (non_zero_term(spins_in)) {
        auto [spins_out, coeff] = term_action(spins_in);
        auto [raw_idx_out, sym, norm_out] = // raw_idx_out = idx_out + 1
            basis_out.representative_data(spins_out);
        if (XDIAG_LIKELY(raw_idx_out)) { // raw_idx_out == 0 means 0-norm
          double norm_in = basis_in.norm(idx_in);
          coeff_t val = coeff * characters(sym) * norm_out / norm_in;
          XDIAG_FILL(idx_in, raw_idx_out - 1, val);
        }
      }
    }
  }
#else
  for (int64_t idx_in = 0; idx_in < basis_in.size(); ++idx_in) {
    bit_t spins_in = basis_in[idx_in];
    if (non_zero_term(spins_in)) {
      auto [spins_out, coeff] = term_action(spins_in);
      auto [raw_idx_out, sym, norm_out] =
          basis_out.representative_data(spins_out);
      if (XDIAG_LIKELY(raw_idx_out)) {
        double norm_in = basis_in.norm(idx_in);
        coeff_t val = coeff * characters(sym) * norm_out / norm_in;
        fill(idx_in, raw_idx_out - 1, val);
      }
    }
  }
#endif
}

} // namespace xdiag::matrices
