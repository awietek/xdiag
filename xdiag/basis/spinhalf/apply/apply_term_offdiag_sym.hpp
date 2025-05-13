// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/common.hpp>
#ifdef _OPENMP
#include <xdiag/parallel/omp/omp_utils.hpp>
#endif

namespace xdiag::basis::spinhalf {

template <typename bit_t, typename coeff_t, class basis_t,
          class non_zero_term_f, class term_action_f, class fill_f>
void apply_term_offdiag_sym_to_spins(bit_t spins_in, int64_t idx_in,
                                     arma::Col<coeff_t> const &characters,
                                     basis_t const &basis_in,
                                     basis_t const &basis_out,
                                     non_zero_term_f non_zero_term,
                                     term_action_f term_action, fill_f fill) {

  if (non_zero_term(spins_in)) {
    auto [spins_out, coeff] = term_action(spins_in);
    auto [idx_out, sym] = basis_out.index_sym(spins_out);

    if (idx_out != invalid_index) {
      double norm_out = basis_out.norm(idx_out);
      double norm_in = basis_in.norm(idx_in);
      coeff_t bloch = characters(sym);
      coeff_t val = coeff * bloch * norm_out / norm_in;
      fill(idx_in, idx_out, val);
    }
  }
}

template <typename coeff_t, class basis_t, class non_zero_term_f,
          class term_action_f, class fill_f>
void apply_term_offdiag_sym(basis_t const &basis_in, basis_t const &basis_out,
                            non_zero_term_f non_zero_term,
                            term_action_f term_action, fill_f fill) {
  using bit_t = typename basis_t::bit_t;

  Representation irrep_out = basis_out.irrep();
  auto characters = irrep_out.characters().as<arma::Col<coeff_t>>();

#ifdef _OPENMP
  int64_t size = basis_in.size();

#pragma omp parallel for schedule(guided)
  for (int64_t idx_in = 0; idx_in < size; ++idx_in) {
    bit_t spins_in = basis_in.state(idx_in);
    apply_term_offdiag_sym_to_spins(spins_in, idx_in, characters, basis_in,
                                    basis_out, non_zero_term, term_action,
                                    fill);
  }

#else
  // Run through all states and apply term
  int64_t idx_in = 0;
  for (auto spins_in : basis_in) {
    apply_term_offdiag_sym_to_spins(spins_in, idx_in, characters, basis_in,
                                    basis_out, non_zero_term, term_action,
                                    fill);
    ++idx_in;
  }
#endif
}

} // namespace xdiag::basis::spinhalf
