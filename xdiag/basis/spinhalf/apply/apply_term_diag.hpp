#pragma once

#include <xdiag/common.hpp>

#ifdef _OPENMP
#include <xdiag/parallel/omp/omp_utils.hpp>
#endif

namespace xdiag::basis::spinhalf {

template <typename bit_t, typename coeff_t, class term_coeff_f, class fill_f>
void apply_term_diag_to_spins(bit_t spins, int64_t idx, term_coeff_f term_coeff,
                              fill_f fill) {
  coeff_t coeff = term_coeff(spins);
  fill(idx, idx, coeff);
}

template <typename coeff_t, class basis_t, class term_coeff_f, class fill_f>
void apply_term_diag(basis_t const &basis, term_coeff_f term_coeff,
                     fill_f fill) {
  using bit_t = typename basis_t::bit_t;

#ifdef _OPENMP
  int64_t size = basis.size();

#pragma omp parallel for schedule(guided)
  for (int64_t idx = 0; idx < size; ++idx) {
    bit_t spins = basis.state(idx);
    apply_term_diag_to_spins<bit_t, coeff_t>(spins, idx, term_coeff, fill);
  }

#else
  int64_t idx = 0;
  for (auto spins : basis) {
    apply_term_diag_to_spins<bit_t, coeff_t>(spins, idx, term_coeff, fill);
    ++idx;
  }
#endif
}

} // namespace xdiag::basis::spinhalf
