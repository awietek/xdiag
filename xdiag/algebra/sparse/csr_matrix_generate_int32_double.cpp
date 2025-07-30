// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "csr_matrix_generate.cpp"

namespace xdiag::algebra {
  template CSRMatrix<int32_t, double>
  csr_matrix_generate<int32_t, double, Spinhalf>(
      OpSum const &, Spinhalf const &, Spinhalf const &, int32_t, bool);

  template CSRMatrix<int32_t, double> csr_matrix_generate<int32_t, double, tJ>(
      OpSum const &, tJ const &, tJ const &, int32_t, bool);

  template CSRMatrix<int32_t, double>
  csr_matrix_generate<int32_t, double, Electron>(
      OpSum const &, Electron const &, Electron const &, int32_t, bool);

} // namespace xdiag::algebra
