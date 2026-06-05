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
#include <xdiag/basis/basis_sublattice.hpp>
#include <xdiag/basis/basis_symmetric.hpp>
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

// Shared implementation for any symmetric basis providing representative_data,
// inv_norm, irrep, size, begin/end (random-access), and operator[].
namespace detail {
template <typename basis_t, typename non_zero_term_f, typename term_action_f,
          typename fill_f>
void term_offdiag_sym(basis_t const &basis_in, basis_t const &basis_out,
                      non_zero_term_f non_zero_term, term_action_f term_action,
                      fill_f fill) {
  using namespace symmetries;
  using bit_t = typename basis_t::bit_t;
  using coeff_t =
      typename std::invoke_result_t<decltype(term_action), bit_t>::second_type;

  auto characters = basis_out.characters().template as<arma::Col<coeff_t>>();

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
        auto [raw_idx_out, sym, norm_out] = // raw_idx_out = idx_out + 1
            basis_out.representative_data(spins_out);
        if (XDIAG_LIKELY(raw_idx_out)) { // raw_idx_out == 0 means 0-norm
          double inv_norm_in = basis_in.inv_norm(idx_in);
          coeff_t val = coeff * characters(sym) * norm_out * inv_norm_in;
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
        double inv_norm_in = basis_in.inv_norm(idx_in);
        coeff_t val = coeff * characters(sym) * norm_out * inv_norm_in;
        fill(idx_in, raw_idx_out - 1, val);
      }
    }
  }
#endif
}
} // namespace detail

template <typename enumeration_t, typename non_zero_term_f,
          typename term_action_f, typename fill_f>
void term_offdiag(basis::BasisSymmetric<enumeration_t> const &basis_in,
                  basis::BasisSymmetric<enumeration_t> const &basis_out,
                  non_zero_term_f non_zero_term, term_action_f term_action,
                  fill_f fill) {
  detail::term_offdiag_sym(basis_in, basis_out, non_zero_term, term_action,
                           fill);
}

template <typename bit_t, int n_sublat, typename non_zero_term_f,
          typename term_action_f, typename fill_f>
void term_offdiag(basis::BasisSublattice<bit_t, n_sublat> const &basis_in,
                  basis::BasisSublattice<bit_t, n_sublat> const &basis_out,
                  non_zero_term_f non_zero_term, term_action_f term_action,
                  fill_f fill) {
  detail::term_offdiag_sym(basis_in, basis_out, non_zero_term, term_action,
                           fill);
}

} // namespace xdiag::matrices
