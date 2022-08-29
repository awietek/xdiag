#pragma once

#include <hydra/common.h>

#ifdef HYDRA_ENABLE_OPENMP
#include <hydra/parallel/omp/omp_utils.h>
#endif

namespace hydra::terms::spinhalf {

template <typename bit_t, typename coeff_t, class IndexingOut,
          class NonZeroTerm, class TermAction, class Fill>
void apply_term_offdiag_no_sym_to_spins(bit_t spins_in, idx_t idx_in,
                                        IndexingOut &&indexing_out,
                                        NonZeroTerm &&non_zero_term,
                                        TermAction &&term_action, Fill &&fill) {
  if (non_zero_term(spins_in)) {
    auto [spins_out, coeff] = term_action(spins_in);
    auto idx_out = indexing_out.index(spins_out);
    assert(idx_out != invalid_index);
    fill(idx_out, idx_in, coeff);
  }
}

template <typename bit_t, typename coeff_t, class IndexingIn, class IndexingOut,
          class NonZeroTerm, class TermAction, class Fill>
void apply_term_offdiag_no_sym(IndexingIn &&indexing_in,
                               IndexingOut &&indexing_out,
                               NonZeroTerm &&non_zero_term,
                               TermAction &&term_action, Fill &&fill) {
#ifdef HYDRA_ENABLE_OPENMP
  idx_t size = indexing_in.size();

#pragma omp parallel for schedule(guided)
  for (idx_t idx_in = 0; idx_in < size; ++idx_in) {
    bit_t spins_in = indexing_in.state(idx_in);
    apply_term_offdiag_no_sym_to_spins<bit_t, coeff_t>(
        spins_in, idx_in, indexing_out, non_zero_term, term_action, fill);
  }

#else
  for (auto [spins_in, idx_in] : indexing_in) {
    apply_term_offdiag_no_sym_to_spins<bit_t, coeff_t>(
        spins_in, idx_in, indexing_out, non_zero_term, term_action, fill);
  }
#endif
}

} // namespace hydra::terms::spinhalf
