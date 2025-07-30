// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "coo_matrix_generate.cpp"

namespace xdiag::algebra {

template COOMatrix<int32_t, complex>
coo_matrix_generate<int32_t, complex, Spinhalf>(OpSum const &, Spinhalf const &,
                                                Spinhalf const &, int32_t);
template COOMatrix<int32_t, complex>
coo_matrix_generate<int32_t, complex, tJ>(OpSum const &, tJ const &, tJ const &,
                                          int32_t);
template COOMatrix<int32_t, complex>
coo_matrix_generate<int32_t, complex, Electron>(OpSum const &, Electron const &,
                                                Electron const &, int32_t);

} // namespace xdiag::algebra
