#pragma once

#include <hydra/common.h>

#ifdef _OPENMP
#include <hydra/parallel/omp/omp_utils.h>
#endif

namespace hydra::spinhalf {

template <typename bit_t, typename coeff_t, class BasisOut,
          class NonZeroTerm, class TermAction, class Fill>
void apply_term_offdiag_no_sym_to_spins(bit_t spins_in, idx_t idx_in,
                                        BasisOut &&basis_out,
                                        NonZeroTerm &&non_zero_term,
                                        TermAction &&term_action, Fill &&fill) {
  if (non_zero_term(spins_in)) {
    auto [spins_out, coeff] = term_action(spins_in);
    auto idx_out = basis_out.index(spins_out);
    assert(idx_out != invalid_index);
    fill(idx_out, idx_in, coeff);
  }
}

template <typename bit_t, typename coeff_t, class BasisIn, class BasisOut,
          class NonZeroTerm, class TermAction, class Fill>
void apply_term_offdiag_no_sym(BasisIn &&basis_in,
                               BasisOut &&basis_out,
                               NonZeroTerm &&non_zero_term,
                               TermAction &&term_action, Fill &&fill) {
#ifdef _OPENMP
  idx_t size = basis_in.size();

#pragma omp parallel for schedule(guided)
  for (idx_t idx_in = 0; idx_in < size; ++idx_in) {
    bit_t spins_in = basis_in.state(idx_in);
    apply_term_offdiag_no_sym_to_spins<bit_t, coeff_t>(
        spins_in, idx_in, basis_out, non_zero_term, term_action, fill);
  }

#else
  for (auto [spins_in, idx_in] : basis_in) {
    apply_term_offdiag_no_sym_to_spins<bit_t, coeff_t>(
        spins_in, idx_in, basis_out, non_zero_term, term_action, fill);
  }
#endif
}

} // namespace hydra::spinhalf
