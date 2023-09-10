#pragma once

#include <hydra/common.h>

#ifdef _OPENMP
#include <hydra/parallel/omp/omp_utils.h>
#endif

namespace hydra::spinhalf {
template <typename bit_t, typename coeff_t, class TermCoeff, class Fill>
void apply_term_diag_to_spins(bit_t spins, int64_t idx, TermCoeff &&term_coeff,
                              Fill &&fill) {
  coeff_t coeff = term_coeff(spins);
  fill(idx, idx, coeff);
}

template <typename bit_t, typename coeff_t, class Basis, class TermCoeff,
          class Fill>
void apply_term_diag(Basis &&basis, TermCoeff &&term_coeff, Fill &&fill) {

#ifdef _OPENMP
  int64_t size = basis.size();

#pragma omp parallel for schedule(guided)
  for (int64_t idx = 0; idx < size; ++idx) {
    bit_t spins = basis.state(idx);
    apply_term_diag_to_spins<bit_t, coeff_t>(spins, idx, term_coeff, fill);
  }

#else
  for (auto [spins, idx] : basis) {
    apply_term_diag_to_spins<bit_t, coeff_t>(spins, idx, term_coeff, fill);
  }
#endif
}

} // namespace hydra::spinhalf
