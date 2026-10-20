// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once
#include <cstdint>

#include <xdiag/operators/coeff.hpp>
#include <xdiag/utils/error.hpp>

namespace xdiag::kernels {

template <typename coeff_t, class basis_t, class fill_f>
void term_identity(Coeff const &coeff, basis_t const &basis, fill_f fill) try {
  coeff_t s = coeff.scalar().as<coeff_t>();
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int64_t i = 0; i < basis.size(); ++i) {
    fill(i, i, s);
  }
}
XDIAG_CATCH

} // namespace xdiag::kernels
