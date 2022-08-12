#pragma once

#include <hydra/common.h>

namespace hydra::terms::spinhalf {

template <typename bit_t, typename coeff_t, class Indexing, class TermCoeff,
          class Fill>
void apply_term_diag(Indexing &&indexing, TermCoeff &&term_coeff, Fill &&fill) {

  for (auto [spins, idx] : indexing) {
    coeff_t coeff = term_coeff(spins);
    fill(idx, idx, coeff);
  }
}

} // namespace hydra::terms::spinhalf
