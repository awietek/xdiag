#pragma once

#include <xdiag/common.hpp>
#ifdef _OPENMP
#include <xdiag/parallel/omp/omp_utils.hpp>
#endif

namespace xdiag::basis::spinhalf {

template <typename bit_t, typename coeff_t, class basis_t,
          class non_zero_term_f, class term_action_f, class fill_f>
void apply_term_offdiag_no_sym_to_spins(bit_t spins_in, int64_t idx_in,
                                        basis_t const &basis_out,
                                        non_zero_term_f non_zero_term,
                                        term_action_f term_action,
                                        fill_f fill) {
  if (non_zero_term(spins_in)) {
    auto [spins_out, coeff] = term_action(spins_in);
    auto idx_out = basis_out.index(spins_out);
    assert(idx_out != invalid_index);
    fill(idx_in, idx_out, coeff);
    // Log("SH in: {} out: {} coeff: {}", idx_in, idx_out, coeff);
  }
}

template <typename bit_t, typename coeff_t, class basis_t,
          class non_zero_term_f, class term_action_f, class fill_f>
void apply_term_offdiag_no_sym(basis_t const &basis_in,
                               basis_t const &basis_out,
                               non_zero_term_f non_zero_term,
                               term_action_f term_action, fill_f fill) {
#ifdef _OPENMP
  int64_t size = basis_in.size();

#pragma omp parallel for schedule(guided)
  for (int64_t idx_in = 0; idx_in < size; ++idx_in) {
    bit_t spins_in = basis_in.state(idx_in);
    apply_term_offdiag_no_sym_to_spins<bit_t, coeff_t>(
        spins_in, idx_in, basis_out, non_zero_term, term_action, fill);
  }

#else
  int64_t idx_in = 0;
  for (auto spins_in : basis_in) {
    apply_term_offdiag_no_sym_to_spins<bit_t, coeff_t>(
        spins_in, idx_in, basis_out, non_zero_term, term_action, fill);
    ++idx_in;
  }
#endif
}

} // namespace xdiag::basis::spinhalf
