// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include <xdiag/algebra/sparse/csr_matrix_nnz.cpp>

namespace xdiag::algebra {
template std::vector<int32_t> csr_matrix_nnz<int32_t, double>(OpSum const &,
                                                              Spinhalf const &,
                                                              Spinhalf const &,
                                                              bool);
template std::vector<int32_t>
csr_matrix_nnz<int32_t, double>(OpSum const &, tJ const &, tJ const &, bool);
template std::vector<int32_t> csr_matrix_nnz<int32_t, double>(OpSum const &,
                                                              Electron const &,
                                                              Electron const &,
                                                              bool);
} // namespace xdiag::algebra
