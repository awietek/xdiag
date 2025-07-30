// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "csr_matrix_generate.cpp"

namespace xdiag::algebra {

template CSRMatrix<int64_t, double>
csr_matrix_generate<int64_t, double, Spinhalf>(OpSum const &, Spinhalf const &,
                                               Spinhalf const &, int64_t, bool);
template CSRMatrix<int64_t, double>
csr_matrix_generate<int64_t, double, tJ>(OpSum const &, tJ const &, tJ const &,
                                         int64_t, bool);
template CSRMatrix<int64_t, double>
csr_matrix_generate<int64_t, double, Electron>(OpSum const &, Electron const &,
                                               Electron const &, int64_t, bool);
} // namespace xdiag::algebra
