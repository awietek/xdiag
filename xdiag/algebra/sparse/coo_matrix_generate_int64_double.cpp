// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "coo_matrix_generate.cpp"

namespace xdiag::algebra {

template COOMatrix<int64_t, double>
coo_matrix_generate<int64_t, double, Spinhalf>(OpSum const &, Spinhalf const &,
                                               Spinhalf const &, int64_t);
template COOMatrix<int64_t, double>
coo_matrix_generate<int64_t, double, tJ>(OpSum const &, tJ const &, tJ const &,
                                         int64_t);
template COOMatrix<int64_t, double>
coo_matrix_generate<int64_t, double, Electron>(OpSum const &, Electron const &,
                                               Electron const &, int64_t);

} // namespace xdiag::algebra
