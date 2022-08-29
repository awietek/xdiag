#pragma once

#include <hydra/common.h>

#ifdef HYDRA_ENABLE_OPENMP
#include <hydra/parallel/omp/omp_utils.h>
#endif

namespace hydra::terms::spinhalf {
template <typename bit_t, typename coeff_t, class TermCoeff, class Fill>
void apply_term_diag_to_spins(bit_t spins, idx_t idx, TermCoeff &&term_coeff,
                              Fill &&fill) {
  coeff_t coeff = term_coeff(spins);
  fill(idx, idx, coeff);
}

template <typename bit_t, typename coeff_t, class Indexing, class TermCoeff,
          class Fill>
void apply_term_diag(Indexing &&indexing, TermCoeff &&term_coeff, Fill &&fill) {

#ifdef HYDRA_ENABLE_OPENMP
    idx_t size = indexing.size();

#pragma omp parallel for schedule(static)
    for (idx_t idx = 0; idx < size; ++idx) {
      bit_t spins = indexing.state(idx);
      apply_term_diag_to_spins<bit_t, coeff_t>(spins, idx, term_coeff, fill);
    }

#else
  for (auto [spins, idx] : indexing) {
    apply_term_diag_to_spins<bit_t, coeff_t>(spins, idx, term_coeff, fill);
  }
#endif
}

} // namespace hydra::terms::spinhalf
