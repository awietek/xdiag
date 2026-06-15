// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <type_traits>
#ifdef _OPENMP
#include <omp.h>
#endif

#include <xdiag/basis/basis_onthefly.hpp>
#include <xdiag/basis/basis_symmetric.hpp>
#include <xdiag/matrices/fill_functions.hpp>
#include <xdiag/matrices/terms/term_offdiag.hpp>
#include <xdiag/utils/likely.hpp>
#include <xdiag/utils/thread_range.hpp>

// Fermionic counterpart of term_offdiag (terms/term_offdiag.hpp). For a
// non-symmetric basis it is identical to term_offdiag: the fermionic sign of
// the operator action is already carried in the term_action coefficient and no
// symmetry permutation is involved. For a symmetric basis it additionally
// accounts for the fermi sign of the permutation that maps the output state to
// its representative; that sign is precomputed per state in the representative
// table, so here it is just a flip of the matrix element when set.

namespace xdiag::matrices {

template <typename enumeration_t, typename non_zero_term_f,
          typename term_action_f, typename fill_f>
void term_offdiag_fermionic(basis::BasisOnTheFly<enumeration_t> const &basis_in,
                            basis::BasisOnTheFly<enumeration_t> const &basis_out,
                            non_zero_term_f non_zero_term,
                            term_action_f term_action, fill_f fill) {
  term_offdiag(basis_in, basis_out, non_zero_term, term_action, fill);
}

template <typename enumeration_t, typename non_zero_term_f,
          typename term_action_f, typename fill_f>
void term_offdiag_fermionic(
    basis::BasisSymmetric<enumeration_t> const &basis_in,
    basis::BasisSymmetric<enumeration_t> const &basis_out,
    non_zero_term_f non_zero_term, term_action_f term_action, fill_f fill) {
  using bit_t = typename enumeration_t::bit_t;
  using coeff_t =
      typename std::invoke_result_t<decltype(term_action), bit_t>::second_type;

  auto characters = basis_out.characters().template as<arma::Col<coeff_t>>();

#ifdef _OPENMP
#pragma omp parallel
  {
    int num_thread = omp_get_thread_num();
    auto [begin, end, idx_in] =
        utils::thread_range(basis_in, num_thread, omp_get_num_threads());
    for (auto it = begin; it != end; ++it, ++idx_in) {
      bit_t spins_in = *it;
      if (non_zero_term(spins_in)) {
        auto [spins_out, coeff] = term_action(spins_in);
        auto [raw_idx_out, sym, norm_out, fermi] =
            basis_out.representative_data_fermi(spins_out);
        if (XDIAG_LIKELY(raw_idx_out)) {
          double inv_norm_in = basis_in.inv_norm(idx_in);
          coeff_t val = coeff * characters(sym) * norm_out * inv_norm_in;
          XDIAG_FILL(idx_in, raw_idx_out - 1, fermi ? -val : val);
        }
      }
    }
  }
#else
  for (int64_t idx_in = 0; idx_in < basis_in.size(); ++idx_in) {
    bit_t spins_in = basis_in[idx_in];
    if (non_zero_term(spins_in)) {
      auto [spins_out, coeff] = term_action(spins_in);
      auto [raw_idx_out, sym, norm_out, fermi] =
          basis_out.representative_data_fermi(spins_out);
      if (XDIAG_LIKELY(raw_idx_out)) {
        double inv_norm_in = basis_in.inv_norm(idx_in);
        coeff_t val = coeff * characters(sym) * norm_out * inv_norm_in;
        fill(idx_in, raw_idx_out - 1, fermi ? -val : val);
      }
    }
  }
#endif
}

} // namespace xdiag::matrices
